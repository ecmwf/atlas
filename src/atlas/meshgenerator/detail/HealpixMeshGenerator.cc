/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <numeric>
#include <vector>

#include "eckit/utils/Hash.h"

#include "atlas/array/Array.h"
#include "atlas/array/ArrayView.h"
#include "atlas/array/IndexView.h"
#include "atlas/field/Field.h"
#include "atlas/grid/Distribution.h"
#include "atlas/grid/Partitioner.h"
#include "atlas/grid/StructuredGrid.h"
#include "atlas/library/config.h"
#include "atlas/mesh/ElementType.h"
#include "atlas/mesh/Elements.h"
#include "atlas/mesh/HybridElements.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/meshgenerator/detail/HealpixMeshGenerator.h"
#include "atlas/meshgenerator/detail/MeshGeneratorFactory.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"
#include "atlas/util/CoordinateEnums.h"
#include "atlas/util/Topology.h"

#define DEBUG_OUTPUT 0
#define DEBUG_OUTPUT_DETAIL 0

using atlas::util::Topology;

namespace atlas {
namespace meshgenerator {

HealpixMeshGenerator::HealpixMeshGenerator(const eckit::Parametrisation& p) {
    std::string mpi_comm = mpi::comm().name();
    p.get("mpi_comm", mpi_comm);
    options.set("mpi_comm",mpi_comm);

    configure_defaults();

    size_t nb_parts;
    if (p.get("nb_parts", nb_parts)) {
        options.set("nb_parts", nb_parts);
    }

    size_t part;
    if (p.get("part", part)) {
        options.set("part", part);
    }

    bool three_dimensional;
    if (p.get("3d", three_dimensional)) {
        options.set("3d", three_dimensional);
    }

    std::string pole_elements{"quads"};
    if (p.get("pole_elements", pole_elements)) {
        if (pole_elements != "pentagons" and pole_elements != "quads") {
            Log::warning() << "Atlas::HealpixMeshGenerator accepts \"pentagons\" or \"quads\" for \"pole_elements\"."
                           << "Defaulting to pole_elements = quads" << std::endl;
        }
    }
    options.set("pole_elements", pole_elements);

    std::string partitioner;
    if (p.get("partitioner", partitioner)) {
        if (not grid::Partitioner::exists(partitioner)) {
            Log::warning() << "Atlas does not have support for partitioner " << partitioner << ". "
                           << "Defaulting to use partitioner EqualRegions" << std::endl;
            partitioner = "equal_regions";
        }
        options.set("partitioner", partitioner);
    }
}

void HealpixMeshGenerator::configure_defaults() {
    std::string mpi_comm;
    options.get("mpi_comm",mpi_comm);
    auto& comm = mpi::comm(mpi_comm);

    // This option sets number of parts the mesh will be split in
    options.set("nb_parts", comm.size());

    // This option sets the part that will be generated
    options.set("part", comm.rank());

    // This option switches between original HEALPix with 1 point at the pole (3d -> true)
    // or HEALPix with 8 or 4 points at the pole (3d -> false)
    options.set("3d", false);

    // This options switches between pentagons and quads as the pole elements for (3d -> false)
    options.set("pole_elements", "quads");

    // This options sets the default partitioner
    std::string partitioner;
    if (grid::Partitioner::exists("equal_regions") && comm.size() > 1) {
        partitioner = "equal_regions";
    }
    else {
        partitioner = "serial";
    }
    options.set<std::string>("partitioner", partitioner);
}


// match glb_idx of node in (nb_pole_nodes==8 and ==4)-meshes to glb_idx of nodes in (nb_pole_nodes==1)-mesh
gidx_t HealpixMeshGenerator::match_node_idx(const gidx_t& gidx, const int ns) const {
    const gidx_t nb_nodes_orig = 12 * ns * ns;
    if (gidx > nb_nodes_ - 1) {
        // no change in index for the periodic nodes
        return gidx;
    }
    if (nb_pole_nodes_ > 1) {
        if (gidx == nb_pole_nodes_ / 2) {
            return 0;
        }
        if (gidx == nb_nodes_ - nb_pole_nodes_ / 2) {
            return nb_nodes_orig + 1;
        }
        bool at_north_pole = (gidx < nb_pole_nodes_);
        bool at_south_pole = (gidx < nb_nodes_ and gidx >= nb_nodes_orig + nb_pole_nodes_);
        if (at_north_pole) {
            return nb_nodes_orig + 2 + gidx - (gidx > nb_pole_nodes_ / 2 ? 1 : 0);
        }
        if (at_south_pole) {
            return gidx - nb_pole_nodes_ + nb_pole_nodes_ + 1 - (gidx > nb_nodes_orig + nb_pole_nodes_ * 3 / 2 ? 1 : 0);
        }
        return gidx - nb_pole_nodes_ + 1;
    }
    // no change for 3d healpix mesh with one node per pole, i.e. nb_pole_nodes = 1
    return gidx;
}


// return "global_id - 1"
gidx_t HealpixMeshGenerator::idx_xy_to_x(const int xidx, const int yidx, const int ns) const {
    ATLAS_ASSERT(yidx < 4 * ns + 1 && yidx >= 0);
    ATLAS_ASSERT(xidx >= 0);

    const gidx_t nb_nodes_orig = 12 * ns * ns;
    auto ghostIdx              = [this](int latid) { return this->nb_nodes_ + latid; };
    gidx_t ret;

    if (yidx == 0) {
        ATLAS_ASSERT(xidx <= nb_pole_nodes_ && xidx >= 0);
        ret = (xidx != nb_pole_nodes_ ? xidx : ghostIdx(yidx));
    }
    else if (yidx < ns) {
        ATLAS_ASSERT(xidx < 4 * yidx + 1 && xidx >= 0);
        ret = (xidx != 4 * yidx ? 2 * yidx * (yidx - 1) + nb_pole_nodes_ + xidx : ghostIdx(yidx));
    }
    else if (yidx <= 2 * ns) {
        ATLAS_ASSERT(xidx < 4 * ns + 1 && xidx >= 0);
        ret = (xidx != 4 * ns ? 2 * ns * (ns - 1) + 4 * ns * (yidx - ns) + nb_pole_nodes_ + xidx : ghostIdx(yidx));
    }
    else if (yidx <= 3 * ns) {
        ATLAS_ASSERT(xidx < 4 * ns + 1 && xidx >= 0);
        ret = (xidx != 4 * ns ? 2 * ns * (3 * ns + 1) + 4 * ns * (yidx - 2 * ns - 1) + nb_pole_nodes_ + xidx
                              : ghostIdx(yidx));
    }
    else if (yidx == 3 * ns + 1 && ns > 1) {
        ATLAS_ASSERT(xidx < 4 * (ns - 1) + 1 && xidx >= 0);
        ret = (xidx != 4 * (ns - 1) ? 2 * ns * (5 * ns + 1) + 4 * ns * (yidx - 3 * ns - 1) + nb_pole_nodes_ + xidx
                                    : ghostIdx(yidx));
    }
    else if (yidx < 4 * ns) {
        ATLAS_ASSERT(xidx < 4 * (ns - (yidx - 3 * ns)) + 1 && xidx >= 0);
        ret =
            (xidx != 4 * (ns - (yidx - 3 * ns)) ? 2 * ns * (5 * ns + 1) + 4 * ns * (yidx - 3 * ns - 1) -
                                                      2 * (yidx - 3 * ns) * (yidx - 3 * ns - 1) + nb_pole_nodes_ + xidx
                                                : ghostIdx(yidx));
    }
    else {
        ATLAS_ASSERT(xidx <= nb_pole_nodes_ && xidx >= 0);
        ret = (xidx != nb_pole_nodes_ ? nb_nodes_orig + nb_pole_nodes_ + xidx : ghostIdx(yidx));
    }
    return ret;
}

// return global_id of the node "above" (xidx,yidx) node
gidx_t HealpixMeshGenerator::up_idx(const int xidx, const int yidx, const int ns) const {
    ATLAS_ASSERT(yidx <= 4 * ns && yidx >= 0);

    const gidx_t nb_nodes_orig = 12 * ns * ns;
    auto ghostIdx              = [this](int latid) { return this->nb_nodes_ + latid; };

    int ret;

    if (yidx == 0) {
        ATLAS_ASSERT(xidx < nb_pole_nodes_);
        ret = (xidx != nb_pole_nodes_ - 1 ? xidx + 1 : ghostIdx(0));
    }
    else if (yidx == 1) {
        ATLAS_ASSERT(xidx < 4);
        ret = (nb_pole_nodes_ == 8 ? 2 * xidx : (nb_pole_nodes_ == 4 ? xidx : 0));
    }
    else if (yidx < ns) {
        ATLAS_ASSERT(xidx < 4 * yidx);
        if (xidx != 4 * yidx - 1) {
            ret = 2 * (yidx - 2) * (yidx - 1) + nb_pole_nodes_ + xidx - std::floor(xidx / (double)yidx);
        }
        else {
            ret = ghostIdx(yidx - 1);
        }
    }
    else if (yidx == ns && ns < 3) {
        ATLAS_ASSERT(xidx < 4 * ns);
        if (xidx != 4 * ns - 1) {
            ret = 2 * ns * (ns - 1) + nb_pole_nodes_ - 4 * (ns - 1) + (xidx + 1) / 2;
        }
        else {
            ret = ghostIdx(yidx - 1);
        }
    }
    else if (yidx == ns) {
        ATLAS_ASSERT(xidx < 4 * ns);
        if (xidx != 4 * ns - 1) {
            ret = 2 * (ns - 2) * (ns - 1) + nb_pole_nodes_ + xidx - std::floor(xidx / (double)yidx);
        }
        else {
            ret = ghostIdx(yidx - 1);
        }
    }
    else if (yidx <= 3 * ns) {
        ATLAS_ASSERT(xidx < 4 * ns);
        int staggering = (yidx - ns) % 2;
        if (xidx != 4 * ns - 1 || (xidx == 4 * ns - 1 && staggering)) {
            ret = 2 * ns * (ns - 1) + nb_pole_nodes_ + 4 * ns * (yidx - ns - 1) + xidx + 1 - staggering;
        }
        else {
            ret = ghostIdx(yidx - 1);
        }
    }
    else if (yidx < 4 * ns - 1) {
        int yidxl = 4 * ns - yidx;
        ATLAS_ASSERT(xidx < 4 * yidxl);
        ret = nb_nodes_orig + nb_pole_nodes_ + 1 - 2 * (yidxl + 2) * (yidxl + 1) + xidx +
              std::floor(xidx / (double)yidxl);
    }
    else if (yidx == 4 * ns - 1) {
        ATLAS_ASSERT(xidx < 4);
        ret = nb_nodes_orig + nb_pole_nodes_ - 3 - (ns == 1 ? 4 : 8) + 2 * xidx;
    }
    else {
        ATLAS_ASSERT(xidx < nb_pole_nodes_);
        if (ns == 1) {
            if (xidx != nb_pole_nodes_ - 1) {
                ret = nb_nodes_orig + nb_pole_nodes_ - 4 +
                      (nb_pole_nodes_ != 4 ? (xidx % 2 ? -4 + (xidx + 1) / 2 : xidx / 2) : xidx);
            }
            else {
                ret = (nb_pole_nodes_ == 4 ? nb_nodes_orig + xidx
                                           : (nb_pole_nodes_ == 8 ? ghostIdx(4 * ns - 2) : ghostIdx(4 * ns)));
            }
        }
        else {
            ret = nb_nodes_orig + nb_pole_nodes_ + (xidx % 2 ? xidx - 12 : xidx - 4 - xidx / 2);
        }
    }
    return ret;
}

// return global_id of the node "below" (xidx,yidx) node
gidx_t HealpixMeshGenerator::down_idx(const int xidx, const int yidx, const int ns) const {
    ATLAS_ASSERT(yidx <= 4 * ns);

    const gidx_t nb_nodes_orig = 12 * ns * ns;
    auto ghostIdx              = [this](int latid) { return this->nb_nodes_ + latid; };

    int ret;

    if (yidx == 0) {
        ATLAS_ASSERT(xidx < nb_pole_nodes_);
        if (ns == 1) {
            if (xidx != nb_pole_nodes_ - 1) {
                if (nb_pole_nodes_ == 8) {
                    ret = 8 + (xidx % 2 ? 4 + (xidx + 1) / 2 : xidx / 2);
                }
                else {
                    ret = 1 + xidx;
                }
            }
            else {
                ret = ghostIdx(2);
            }
        }
        else {
            ret = nb_pole_nodes_ + ((xidx + 1) % 2 ? xidx / 2 : 4 + xidx);
        }
    }
    else if (yidx < ns) {
        ATLAS_ASSERT(xidx < 4 * yidx);
        ret = 2 * yidx * (yidx + 1) + nb_pole_nodes_ + 1 + xidx + std::floor(xidx / (double)yidx);
    }
    else if (yidx == ns && ns == 1) {
        ATLAS_ASSERT(xidx < 4);
        ret = (xidx != 3 ? nb_pole_nodes_ + 5 + xidx : ghostIdx(2));
    }
    else if (yidx == 2 * ns && ns == 1) {
        ATLAS_ASSERT(xidx < 4);
        ret = (nb_pole_nodes_ == 8 ? 16 + xidx : (nb_pole_nodes_ == 4 ? 12 + xidx : 9 + xidx));
    }
    else if (yidx < 3 * ns && ns > 1) {
        ATLAS_ASSERT(xidx < 4 * ns);
        int staggering = (yidx - ns) % 2;
        if (xidx != 4 * ns - 1 || (xidx == 4 * ns - 1 && staggering)) {
            ret = 2 * ns * (ns - 1) + nb_pole_nodes_ + 4 * ns * (yidx - ns + 1) + xidx +
                  (yidx != 3 * ns ? 1 - staggering : 0);
        }
        else {
            ret = ghostIdx(yidx + 1);
        }
    }
    else if (yidx == 4 * ns - 2) {
        ATLAS_ASSERT(xidx < 8);
        if (nb_pole_nodes_ == 8) {
            ret = (xidx != 7 ? nb_nodes_orig + 4 + (xidx + 1) / 2 : ghostIdx(4 * ns - 1));
        }
        else if (nb_pole_nodes_ == 4) {
            ret = (xidx != 7 ? nb_nodes_orig + (xidx + 1) / 2 : ghostIdx(4 * ns - 1));
        }
        else {
            ret = (xidx != 7 ? nb_nodes_orig - 3 + (xidx + 1) / 2 : ghostIdx(4 * ns - 1));
        }
    }
    else if (yidx == 4 * ns - 1) {
        ATLAS_ASSERT(xidx < 4);
        if (nb_pole_nodes_ == 8) {
            ret = nb_nodes_orig + nb_pole_nodes_ + 2 * xidx;
        }
        else if (nb_pole_nodes_ == 4) {
            ret = nb_nodes_orig + nb_pole_nodes_ + xidx;
        }
        else {
            ret = nb_nodes_orig + 1;
        }
    }
    else if (yidx < 4 * ns - 1) {
        int yidxl = yidx - 3 * ns;
        ATLAS_ASSERT(xidx < 4 * (ns - yidxl));
        if (xidx != 4 * (ns - yidxl) - 1) {
            ret = 2 * ns * (5 * ns + 1) + nb_pole_nodes_ + 4 * ns * yidxl - 2 * (yidxl + 1) * yidxl + xidx -
                  std::floor(xidx / (double)(ns - yidxl));
        }
        else {
            ret = ghostIdx(yidx + 1);
        }
    }
    else if (yidx == 4 * ns) {
        ATLAS_ASSERT(xidx < nb_pole_nodes_);
        ret = (xidx != nb_pole_nodes_ - 1 ? nb_nodes_orig + nb_pole_nodes_ + xidx + 1 : ghostIdx(yidx));
    }
    else {
        throw_AssertionFailed("Invalid value of yidx", Here());
    }
    return ret;
}

// return global_id of the node "to the right of" (xidx,yidx) node
gidx_t HealpixMeshGenerator::right_idx(const int xidx, const int yidx, const int ns) const {
    ATLAS_ASSERT(yidx <= 4 * ns);

    const gidx_t nb_nodes_orig = 12 * ns * ns;
    auto ghostIdx              = [this](int latid) { return this->nb_nodes_ + latid; };
    int ret                    = -1;

    if (yidx == 0) {
        ATLAS_ASSERT(xidx < nb_pole_nodes_);
        if (ns == 1) {
            ret = (xidx != nb_pole_nodes_ - 1
                       ? (xidx % 2 ? nb_pole_nodes_ + (xidx + 1) / 2 : nb_pole_nodes_ + 5 + xidx / 2)
                       : ghostIdx(1));
        }
        else {
            ret = (xidx != nb_pole_nodes_ - 1 ? (xidx % 2 ? nb_pole_nodes_ + (xidx + 1) / 2 : nb_pole_nodes_ + 5 + xidx)
                                              : ghostIdx(1));
        }
    }
    else if (yidx == 1) {
        ATLAS_ASSERT(xidx < 4);
        if (nb_pole_nodes_ == 8) {
            ret = 1 + 2 * xidx;
        }
        else if (nb_pole_nodes_ == 4) {
            ret = (xidx != 3 ? xidx + 1 : ghostIdx(0));
        }
        else {
            ret = (xidx != 3 ? xidx + 2 : ghostIdx(1));
        }
    }
    else if (yidx < ns) {
        ATLAS_ASSERT(xidx < 4 * yidx);
        ret = (xidx != 4 * yidx - 1 ? 2 * yidx * (yidx - 1) + nb_pole_nodes_ + 1 + xidx : ghostIdx(yidx));
    }
    else if (yidx == 3 && ns == 1) {
        ATLAS_ASSERT(xidx < 4);
        if (nb_pole_nodes_ == 8) {
            ret = 21 + 2 * xidx;
        }
        else if (nb_pole_nodes_ == 4) {
            ret = (xidx != 3 ? 17 + xidx : ghostIdx(yidx + 1));
        }
        else {
            ret = (xidx != 3 ? 10 + xidx : ghostIdx(yidx));
        }
    }
    else if (yidx <= 3 * ns) {
        ATLAS_ASSERT(xidx < 4 * ns + 1);
        ret = (xidx != 4 * ns - 1 ? 2 * ns * (ns - 1) + 4 * ns * (yidx - ns) + nb_pole_nodes_ + 1 + xidx
                                  : ghostIdx(yidx));
    }
    else if (yidx < 4 * ns - 1 && ns > 1) {
        int yidxl = yidx - 3 * ns;
        ATLAS_ASSERT(xidx < 4 * (ns - yidxl));
        if (xidx != 4 * (ns - yidxl) - 1) {
            ret = 2 * ns * (5 * ns + 1) + 4 * ns * (yidx - 3 * ns - 1) - 2 * (yidx - 3 * ns) * (yidx - 3 * ns - 1) +
                  nb_pole_nodes_ + 1 + xidx;
        }
        else {
            ret = ghostIdx(yidx);
        }
    }
    else if (yidx == 4 * ns - 1) {
        ATLAS_ASSERT(xidx < 4);
        if (nb_pole_nodes_ == 8) {
            ret = nb_nodes_orig + nb_pole_nodes_ + 1 + 2 * xidx;
        }
        else if (nb_pole_nodes_ == 4) {
            ret = (xidx != 3 ? nb_nodes_orig + nb_pole_nodes_ + 1 + xidx : ghostIdx(yidx + 1));
        }
        else {
            ret = (xidx != 3 ? nb_nodes_orig - 2 + xidx : ghostIdx(yidx));
        }
    }
    else if (yidx == 4 * ns) {
        ATLAS_ASSERT(xidx < nb_pole_nodes_);
        if (xidx != nb_pole_nodes_ - 1) {
            ret = (xidx % 2 ? nb_nodes_orig + 4 + (xidx + 1) / 2 : nb_nodes_orig + 4 - (ns == 1 ? 3 : 7) + xidx);
        }
        else {
            ret = ghostIdx(yidx - 1);
        }
    }
    return ret;
}

// return global_id - 1 of the pentagon node "to the right of" (xidx,yidx) node
// pentagon points are only needed for yidx == 1 and yidx == 4 * ns - 1
gidx_t HealpixMeshGenerator::pentagon_right_idx(const int xidx, const int yidx, const int ns) const {
    auto ghostIdx = [this](int latid) { return this->nb_nodes_ + latid; };
    if (yidx == 1) {
        return (xidx != 3 ? nb_pole_nodes_ + xidx + 1 : ghostIdx(1));
    }
    else if (yidx == 4 * ns - 1) {
        return (xidx != 3 ? (12 * ns * ns + xidx + 1) : ghostIdx(yidx));
    }
    else {
        return -2;
    }
}

void HealpixMeshGenerator::generate(const Grid& grid, Mesh& mesh) const {
    ATLAS_ASSERT(HealpixGrid(grid), "Grid could not be cast to a HealpixGrid");
    ATLAS_ASSERT(!mesh.generated());

    const StructuredGrid rg = StructuredGrid(grid);
    if (!rg) {
        throw_Exception("HealpixMeshGenerator can only work with a Healpix grid", Here());
    }

    size_t nb_parts = options.get<size_t>("nb_parts");

    std::string partitioner_type = "equal_regions";
    options.get("partitioner", partitioner_type);

    mpi::push(options.getString("mpi_comm"));
    grid::Partitioner partitioner(partitioner_type, nb_parts);
    grid::Distribution distribution(partitioner.partition(grid));
    mpi::pop();
    generate(grid, distribution, mesh);
}

void HealpixMeshGenerator::hash(eckit::Hash& h) const {
    h.add("HealpixMeshGenerator");
    options.hash(h);
}

void HealpixMeshGenerator::generate(const Grid& grid, const grid::Distribution& distribution, Mesh& mesh) const {
    ATLAS_TRACE();
    Log::debug() << "HealpixMeshGenerator generating mesh from " << grid.name() << std::endl;
    ATLAS_ASSERT(HealpixGrid(grid), "Grid could not be cast to a HealpixGrid");
    ATLAS_ASSERT(!mesh.generated());

    if (grid.size() != static_cast<idx_t>(distribution.size())) {
        std::stringstream msg;
        msg << "Number of points in grid (" << grid.size()
            << ") different from "
               "number of points in grid distribution ("
            << distribution.size() << ")";
        throw_AssertionFailed(msg.str(), Here());
    }

    // clone some grid properties
    setGrid(mesh, grid, distribution);

    generate_mesh(grid, distribution, mesh);
}

void HealpixMeshGenerator::generate_mesh(const StructuredGrid& grid, const grid::Distribution& distribution,
                                         Mesh& mesh) const {
    // This function should do the following:
    // - define nodes with
    //      mesh.nodes().resize(nnodes);
    //      mesh::Nodes& nodes = mesh.nodes();
    //    following properties should be defined:
    //      array::ArrayView<double,2> xy            ( nodes.xy() );
    //      array::ArrayView<gidx_t,1> glb_idx       ( nodes.global_index() );
    //      array::ArrayView<int,   1> part          ( nodes.partition() );
    //      array::ArrayView<int,   1> ghost         ( nodes.ghost() );
    //      array::ArrayView<int,   1> flags         ( nodes.flags() );
    // - define cells (only quadrilaterals for now) with
    //      mesh.cells().add( mesh::ElementType::create("Quadrilateral"), nquads  );
    //      mesh.cells().add( mesh::ElementType::create("Pentagon"), npents  );
    //    further define cells with
    //      array::ArrayView<gidx_t,1> cells_glb_idx( mesh.cells().global_index()
    //      );
    //      array::ArrayView<int,1>    cells_part(    mesh.cells().partition() );
    // - define connectivity with
    //      mesh::HybridElements::Connectivity& node_connectivity =
    //      mesh.cells().node_connectivity();
    //      node_connectivity.set( jcell, cell_nodes );
    //    where cell_nodes is a 4-element or 5-element integer array containing the LOCAL
    //    indices of the nodes
    //
    // The rule do determine if a cell belongs to a proc is the following with some
    // exceptions near poles: cell belongs to a proc if the left corner of the cell belongs to that proc.


    ATLAS_TRACE();

    ATLAS_ASSERT(HealpixGrid(grid));

    const int mypart                = options.get<size_t>("part");
    const int nparts                = options.get<size_t>("nb_parts");
    const bool three_dimensional    = options.get<bool>("3d");
    const std::string pole_elements = options.get<std::string>("pole_elements");
    const int nb_pole_nodes         = (pole_elements == "pentagons") ? 4 : (three_dimensional ? 1 : 8);
    const int ny                    = grid.ny() + 2;
    const int ns                    = (ny - 1) / 4;
    const int nvertices             = 12 * ns * ns + 2 * nb_pole_nodes;


    nb_pole_nodes_ = nb_pole_nodes;
    nb_points_     = 12 * ns * ns + (nb_pole_nodes == 8 ? 8 : 0);
    nb_nodes_      = nvertices;

    int inode;
    auto nb_lat_nodes = [ny, nb_pole_nodes, &grid](int latid) {
        return ((latid == 0) or (latid == ny - 1) ? nb_pole_nodes : grid.nx()[latid - 1]);
    };

    int ii, ix, iy, ii_ghost, ii_glb;
    int iy_min, iy_max;   // a belt (iy_min:iy_max) surrounding the nodes on this processor
    int nnodes_nonghost;  // non-ghost node: belongs to this part

    // vector of local indices: necessary for remote indices of ghost nodes
    std::vector<int> local_idx(nvertices, -1);
    std::vector<int> current_idx(nparts, 0);  // index counter for each proc

    // ANSATZ: requirement on the partitioner
    auto compute_part = [&](int iy, gidx_t ii_glb) -> int {
        // nodes at the pole belong to proc_0 (north) and proc_maxRank (south)
        // a node gets its proc rank from the element for which this node would be its west vertex
        return (iy == 0 ? 0 : (iy == ny - 1 ? nparts - 1 : distribution.partition(ii_glb - nb_pole_nodes)));
    };

#if DEBUG_OUTPUT_DETAIL
    for (iy = 0; iy < ny; iy++) {
        int nx = nb_lat_nodes(iy);
        for (ix = 0; ix < nx; ix++) {
            Log::info() << "iy, ix, glb_idx, up_idx, down_idx, right_idx, pent_right_idx : " << iy << ", " << ix << ", "
                        << idx_xy_to_x(ix, iy, ns) + 1 << ", " << up_idx(ix, iy, ns) + 1 << ", "
                        << down_idx(ix, iy, ns) + 1 << ", " << right_idx(ix, iy, ns) + 1 << ", "
                        << pentagon_right_idx(ix, iy, ns) + 1 << std::endl;
        }
        Log::info() << std::endl;
    }
#endif

    // loop over all points to determine local indices and surrounding rectangle
    ii_glb          = 0;  // global index starting from 0
    iy_min          = ny + 1;
    iy_max          = 0;
    nnodes_nonghost = 0;
    for (iy = 0; iy < ny; iy++) {
        int nx = nb_lat_nodes(iy);
        for (ix = 0; ix < nx; ix++) {
            int proc_id       = compute_part(iy, ii_glb);
            local_idx[ii_glb] = current_idx[proc_id]++;
            if (proc_id == mypart) {
                ++nnodes_nonghost;
                iy_min = std::min(iy_min, iy);
                iy_max = std::max(iy_max, iy);
            }
            ++ii_glb;
        }
    }

#if DEBUG_OUTPUT_DETAIL
    inode = 0;
    Log::info() << "local_idx : " << std::endl;
    for (size_t ilat = 0; ilat < ny; ilat++) {
        for (size_t ilon = 0; ilon < nb_lat_nodes(ilat); ilon++) {
            Log::info() << std::setw(4) << local_idx[inode];
            inode++;
        }
        Log::info() << std::endl;
    }
    inode = 0;
    Log::info() << "global_idx : " << std::endl;
    for (size_t ilat = 0; ilat < ny; ilat++) {
        for (size_t ilon = 0; ilon < nb_lat_nodes(ilat); ilon++) {
            Log::info() << std::setw(4) << inode;
            inode++;
        }
        Log::info() << std::endl;
    }
#endif

    // dimensions of surrounding belt (SB)
    int nnodes_SB = 0;
    if (iy_min <= 2) {
        iy_min = 0;
    }
    else {
        --iy_min;
    }
    if (iy_max >= ny - 3) {
        iy_max = ny - 1;
    }
    else {
        ++iy_max;
    }
    for (int iy = iy_min; iy <= iy_max; iy++) {
        // (east) periodic point adds +1 here
        nnodes_SB += nb_lat_nodes(iy) + 1;
    }

#if DEBUG_OUTPUT
    Log::info() << "[" << mypart << "] : nnodes_SB = " << nnodes_SB << "\n";
#endif

    // partitions and local indices in SB
    std::vector<int> parts_SB(nnodes_SB, -1);
    std::vector<int> local_idx_SB(nnodes_SB, -1);
    std::vector<bool> is_ghost_SB(nnodes_SB, true);

    // starting from index 0, first global node-index for this partition
    int parts_sidx = idx_xy_to_x(0, iy_min, ns);

    ii       = 0;                                  // index inside SB
    ii_ghost = nnodes_SB - (iy_max - iy_min + 1);  // first local ghost idx
    for (iy = iy_min; iy <= iy_max; iy++) {
        int nx = nb_lat_nodes(iy) + 1;
        for (ix = 0; ix < nx; ix++) {
            if (ix != nx - 1) {
                ii_glb           = ii + parts_sidx;
                parts_SB[ii]     = compute_part(iy, ii_glb);
                local_idx_SB[ii] = ii;
                is_ghost_SB[ii]  = !(parts_SB[ii] == mypart);
                ++ii;
            }
            else {
                parts_SB[ii_ghost]     = compute_part(iy, ii_glb);
                local_idx_SB[ii_ghost] = ii_ghost;
                is_ghost_SB[ii_ghost]  = true;
                ++ii_ghost;
            }
        }
    }

#if DEBUG_OUTPUT_DETAIL
    Log::info() << "[" << mypart << "] : "
                << "parts_SB = ";
    for (ii = 0; ii < nnodes_SB; ii++) {
        Log::info() << parts_SB[ii] << ",";
    }
    Log::info() << std::endl;
    Log::info() << "[" << mypart << "] : "
                << "local_idx_SB = ";
    for (ii = 0; ii < nnodes_SB; ii++) {
        Log::info() << local_idx_SB[ii] << ",";
    }
    Log::info() << std::endl;
    Log::info() << "[" << mypart << "] : "
                << "is_ghost_SB = ";
    for (ii = 0; ii < nnodes_SB; ii++) {
        Log::info() << is_ghost_SB[ii] << ",";
    }
    Log::info() << std::endl;
#endif

    // vectors marking nodes that are necessary for this proc's cells
    std::vector<bool> is_node_SB(nnodes_SB, false);

    // determine number of cells and number of nodes
    int nnodes = 0;
    int nquads = 0;
    int npents = 0;
    ii         = 0;

    int glb2loc_ghost_offset = -nnodes_SB + iy_max + nb_nodes_ + 1;
    auto get_local_id        = [this, &parts_sidx, &glb2loc_ghost_offset](gidx_t gidx) {
        return gidx - (gidx < nb_nodes_ ? parts_sidx : glb2loc_ghost_offset);
    };

    for (iy = iy_min; iy <= iy_max; iy++) {
        int nx = nb_lat_nodes(iy);
        for (ix = 0; ix < nx; ix++) {
            if (not is_ghost_SB[ii]) {
                const bool at_pole      = (iy == 0 or iy == ny - 1);
                bool not_duplicate_cell = (at_pole ? (ix % 2) : 1);
                if (at_pole and nb_pole_nodes < 8) {
                    // nodes at the poles do not own any pentagons
                    // if nb_pole_node=1, the node at the poles does not own any quads
                    not_duplicate_cell = false;
                }
                if (not at_pole or not_duplicate_cell) {
                    if (not is_node_SB[ii]) {
                        ++nnodes;
                        is_node_SB[ii] = true;
                    }
                    bool is_pentagon = (pole_elements == "pentagons" and (iy == 1 or iy == ny - 2));
                    if (is_pentagon) {
                        ++npents;
                        // mark the pentagon right node
                        idx_t iil = get_local_id(pentagon_right_idx(ix, iy, ns));
                        if (not is_node_SB[iil]) {
                            ++nnodes;
                            is_node_SB[iil] = true;
                        }
                    }
                    else {
                        ++nquads;
                    }
                    // mark upper corner
                    idx_t iil = get_local_id(up_idx(ix, iy, ns));
                    if (not is_node_SB[iil]) {
                        ++nnodes;
                        is_node_SB[iil] = true;
                    }
                    // mark lower corner
                    iil = get_local_id(down_idx(ix, iy, ns));
                    if (not is_node_SB[iil]) {
                        ++nnodes;
                        is_node_SB[iil] = true;
                    }
                    // mark right corner
                    iil = get_local_id(right_idx(ix, iy, ns));
                    if (not is_node_SB[iil]) {
                        ++nnodes;
                        is_node_SB[iil] = true;
                    }
                }
            }
            ++ii;
        }
    }
    int ncells = nquads + npents;
    ATLAS_ASSERT(ncells > 0);

#if DEBUG_OUTPUT
    Log::info() << "[" << mypart << "] : "
                << "nnodes = " << nnodes << ", nquads = " << nquads << ", npents = " << npents
                << ", parts_sidx = " << parts_sidx << std::endl;
    Log::info() << "[" << mypart << "] : "
                << "iy_min = " << iy_min << ", iy_max = " << iy_max << std::endl;
#endif
#if DEBUG_OUTPUT_DETAIL
    Log::info() << "[" << mypart << "] : "
                << "is_node_SB = ";
    for (int ii = 0; ii < nnodes_SB; ii++) {
        Log::info() << is_node_SB[ii] << ",";
    }
    Log::info() << std::endl;
#endif

    // define nodes and associated properties
    mesh.nodes().resize(nnodes);
    mesh::Nodes& nodes = mesh.nodes();
    auto xy            = array::make_view<double, 2>(nodes.xy());
    auto lonlat        = array::make_view<double, 2>(nodes.lonlat());
    auto glb_idx       = array::make_view<gidx_t, 1>(nodes.global_index());
    auto remote_idx    = array::make_indexview<idx_t, 1>(nodes.remote_index());
    auto part          = array::make_view<int, 1>(nodes.partition());
    auto ghost         = array::make_view<int, 1>(nodes.ghost());
    auto halo          = array::make_view<int, 1>(nodes.halo());
    auto flags         = array::make_view<int, 1>(nodes.flags());

    int quad_begin;
    int pent_begin;

    // define cells and associated properties
    mesh.cells().add(mesh::ElementType::create("Quadrilateral"), nquads);
    quad_begin = mesh.cells().elements(0).begin();
    pent_begin = quad_begin + mesh.cells().elements(0).size();
    if (pole_elements == "pentagons") {
        mesh.cells().add(mesh::ElementType::create("Pentagon"), npents);
        pent_begin = mesh.cells().elements(1).begin();
    }
    auto cells_part         = array::make_view<int, 1>(mesh.cells().partition());
    auto cells_glb_idx      = array::make_view<gidx_t, 1>(mesh.cells().global_index());
    auto& node_connectivity = mesh.cells().node_connectivity();

    idx_t cell_nodes[5];
    int jquadcell = quad_begin;
    int jpentcell = pent_begin;

    int inode_nonghost, inode_ghost;

    // loop over nodes and set properties
    inode_nonghost = 0;
    inode_ghost    = nnodes_nonghost;  // ghost nodes start counting after nonghost nodes
    ii             = 0;                // index inside SB
    for (iy = iy_min; iy <= iy_max; iy++) {
        int nx = nb_lat_nodes(iy) + 1;
        for (ix = 0; ix < nx; ix++) {
            int iil = idx_xy_to_x(ix, iy, ns);
            iil -= (iil < nb_nodes_ ? parts_sidx : -nnodes_SB + iy_max + nb_nodes_ + 1);
            if (is_node_SB[iil]) {
                // set node counter
                if (is_ghost_SB[iil]) {
                    inode = inode_ghost++;
                }
                else {
                    inode = inode_nonghost++;
                }

                // flags
                Topology::reset(flags(inode));

                glb_idx(inode) = 1 + match_node_idx(idx_xy_to_x(ix, iy, ns), ns);

                // grid coordinates
                double _xy[2];
                double xy1[2], xy2[2];
                if (iy == 0) {
                    _xy[0] = (nb_pole_nodes == 8 ? 45. * ix : (nb_pole_nodes == 4 ? 90. * ix : 180.));
                    _xy[1] = 90.;
                    Topology::set(flags(inode), Topology::BC);
                }
                else if (iy == ny - 1) {
                    _xy[0] = (nb_pole_nodes == 8 ? 45. * ix : (nb_pole_nodes == 4 ? 90. * ix : 180.));
                    _xy[1] = -90.;
                    Topology::set(flags(inode), Topology::BC);
                }
                else if (ix == nx - 1) {
                    grid.xy(ix - 1, iy - 1, xy1);
                    grid.xy(ix - 2, iy - 1, xy2);
                    _xy[0] = 1.5 * xy1[0] - 0.5 * xy2[0];
                    _xy[1] = xy1[1];
                    Topology::set(flags(inode), Topology::BC);
                }
                else if (ix == 0) {
                    grid.xy(ix + 1, iy - 1, xy1);
                    grid.xy(ix, iy - 1, xy2);
                    _xy[0] = 1.5 * xy2[0] - 0.5 * xy1[0];
                    _xy[1] = xy1[1];
                    Topology::set(flags(inode), Topology::BC);
                }
                else {
                    grid.xy(ix, iy - 1, xy1);
                    grid.xy(ix - 1, iy - 1, xy2);
                    _xy[0] = 0.5 * (xy1[0] + xy2[0]);
                    _xy[1] = xy1[1];
                }

                if (Topology::check(flags(inode), Topology::BC)) {
                    if (iy == 0) {
                        Topology::set(flags(inode), Topology::NORTH);
                    }
                    else if (iy == ny - 1) {
                        Topology::set(flags(inode), Topology::SOUTH);
                    }
                    if (ix == 0) {
                        Topology::set(flags(inode), Topology::WEST);
                    }
                    else if (ix == nx - 1) {
                        Topology::set(flags(inode), Topology::EAST | Topology::GHOST);
                        ATLAS_ASSERT(is_ghost_SB[iil]);
                    }
                }

                xy(inode, LON) = _xy[LON];
                xy(inode, LAT) = _xy[LAT];

                // geographic coordinates by using projection
                grid.projection().xy2lonlat(_xy);
                lonlat(inode, LON) = _xy[LON];
                lonlat(inode, LAT) = _xy[LAT];

                part(inode)  = parts_SB[iil];
                ghost(inode) = is_ghost_SB[iil];
                halo(inode)  = 0;

                if (ghost(inode)) {
                    Topology::set(flags(inode), Topology::GHOST);
                    remote_idx(inode) = local_idx_SB[iil];
                }
                else {
                    remote_idx(inode) = -1;
                }
                if (Topology::check(flags(inode), Topology::BC | Topology::EAST)) {
                    part(inode) = mypart;  // To be fixed later
                }
                local_idx_SB[iil] = inode;

#if DEBUG_OUTPUT_DETAIL
                Log::info() << "[" << mypart << "] : "
                            << "New node \tinode=" << inode << "; iil= " << iil << "; ix=" << ix << "; iy=" << iy
                            << "; glon=" << lonlat(inode, 0) << "; glat=" << lonlat(inode, 1)
                            << "; glb_idx=" << glb_idx(inode) << "; loc_idx=" << local_idx_SB[iil] << std::endl;
#endif
            }
            ii += (ix != nx - 1 ? 1 : 0);
        }
    }

    auto get_local_idx_SB = [&local_idx_SB, &get_local_id](gidx_t gidx) { return local_idx_SB[get_local_id(gidx)]; };

    ii               = 0;  // index inside SB (surrounding belt)
    int jcell_offset = 0;  // global index offset due to extra points at the north pole
    gidx_t jcell     = 0;  // global cell counter
    gidx_t points_in_partition = 0;
    for (iy = iy_min; iy <= iy_max; iy++) {
        int nx = nb_lat_nodes(iy) + 1;
        points_in_partition += (nx - 1);
        for (ix = 0; ix < nx; ix++) {
            const bool at_pole     = (iy == 0 or iy == ny - 1);
            int not_duplicate_cell = (at_pole ? ix % 2 : 1);
            if (at_pole and (not not_duplicate_cell or nb_pole_nodes < 8)) {
                // if nb_pole_nodes = 4 : the pole nodes do not own any pentagons
                // if nb_pole_nodes = 1 : the pole nodes do not own any quads
                continue;
            }
            bool pentagon = (iy == 1 or iy == 4 * ns - 1) and pole_elements == "pentagons";
            int iil       = idx_xy_to_x(ix, iy, ns);
            iil -= (iil < nb_nodes_ ? parts_sidx : -nnodes_SB + iy_max + nb_nodes_ + 1);
            if (not is_ghost_SB[iil]) {
                jcell++;  // a cell will be added
                // define cell vertices (in local indices) in cell_nodes
                int node_id           = 0;
                cell_nodes[node_id++] = get_local_idx_SB(idx_xy_to_x(ix, iy, ns));
                cell_nodes[node_id++] = get_local_idx_SB(down_idx(ix, iy, ns));
                bool south_hemisphere = (iy > 2 * ns);
                if (pentagon && not south_hemisphere) {
                    cell_nodes[node_id++] = get_local_idx_SB(pentagon_right_idx(ix, iy, ns));
                }
                cell_nodes[node_id++] = get_local_idx_SB(right_idx(ix, iy, ns));
                if (pentagon && south_hemisphere) {
                    // in the south hemisphere the pentagon point comes in a different clock-wise ordering
                    cell_nodes[node_id++] = get_local_idx_SB(pentagon_right_idx(ix, iy, ns));
                }
                cell_nodes[node_id++] = get_local_idx_SB(up_idx(ix, iy, ns));

                // match global cell indexing for the three healpix versions
                if (nb_pole_nodes == 1) {
                    cells_glb_idx(jquadcell) = parts_sidx + points_in_partition - nx + ix + 1;
                }
                else if (nb_pole_nodes == 8) {
                    if (iy == 0) {
                        cells_glb_idx(jquadcell) = 12 * ns * ns + 1 + ix / 2;
                        jcell_offset++;
                    }
                    else if (iy == ny - 1) {
                        cells_glb_idx(jquadcell) = 12 * ns * ns + 5 + ix / 2;
                        jcell_offset++;
                    }
                    else {
                        cells_glb_idx(jquadcell) = parts_sidx + iil - 3 - (mypart != 0 ? 4 : jcell_offset);
                    }
                }
                else if (nb_pole_nodes == 4) {
                    if (iy == 0 or iy == ny - 1) {
                        continue;
                    }
                    if (pentagon) {
                        cells_glb_idx(jpentcell) = (mypart != 0 ? 12 * ns * ns - 3 + ix : jcell);
                        jcell_offset++;
                    }
                    else {
                        cells_glb_idx(jquadcell) = parts_sidx + points_in_partition - nx - 6 + ix + (mypart != 0 ? 4 : jcell_offset);
                    }
                }
#if DEBUG_OUTPUT_DETAIL
                std::cout << "[" << mypart << "] : ";
                if (pentagon) {
                    std::cout << "New pent: loc-idx " << jpentcell << ", glb-idx " << cells_glb_idx(jpentcell) << ": ";
                }
                else {
                    std::cout << "New quad: loc-idx " << jquadcell << ", glb-idx " << cells_glb_idx(jquadcell) << ": ";
                }
                std::cout << glb_idx(cell_nodes[0]) << "," << glb_idx(cell_nodes[1]) << "," << glb_idx(cell_nodes[2])
                          << "," << glb_idx(cell_nodes[3]);
                if (pentagon) {
                    std::cout << "," << glb_idx(cell_nodes[4]);
                }
                std::cout << std::endl;
#endif
                // add cell to the node connectivity table
                if (pentagon) {
                    cells_part(jpentcell) = mypart;
                    node_connectivity.set(jpentcell, cell_nodes);
                    ++jpentcell;
                }
                else {
                    cells_part(jquadcell) = mypart;
                    node_connectivity.set(jquadcell, cell_nodes);
                    ++jquadcell;
                }
            }
            ii += (ix != nx - 1 ? 1 : 0);
        }
    }

#if DEBUG_OUTPUT_DETAIL
    // list nodes
    Log::info() << "Listing nodes ...";
    for (inode = 0; inode < nnodes; inode++) {
        std::cout << "[" << mypart << "] : "
                  << " node " << inode << ": ghost = " << ghost(inode) << ", glb_idx = " << glb_idx(inode)
                  << ", part = " << part(inode) << ", lon = " << lonlat(inode, 0) << ", lat = " << lonlat(inode, 1)
                  << ", remote_idx = " << remote_idx(inode) << std::endl;
    }

    for (gidx_t jcell = quad_begin; jcell < nquads; jcell++) {
        std::cout << "[" << mypart << "] : "
                  << " cell " << jcell << ", glb-idx " << cells_glb_idx(jcell) << ": "
                  << glb_idx(node_connectivity(jcell, 0)) << "," << glb_idx(node_connectivity(jcell, 1)) << ","
                  << glb_idx(node_connectivity(jcell, 2)) << "," << glb_idx(node_connectivity(jcell, 3)) << std::endl;
    }
    for (gidx_t jcell = pent_begin; jcell < nquads; jcell++) {
        std::cout << "[" << mypart << "] : "
                  << " cell " << jcell << ", glb-idx " << cells_glb_idx(jcell) << ": "
                  << glb_idx(node_connectivity(jcell, 0)) << "," << glb_idx(node_connectivity(jcell, 1)) << ","
                  << glb_idx(node_connectivity(jcell, 2)) << "," << glb_idx(node_connectivity(jcell, 3)) << std::endl;
    }
#endif

    mesh.metadata().set("nb_parts",options.getInt("nb_parts"));
    mesh.metadata().set("part",options.getInt("part"));
    mesh.metadata().set("mpi_comm",options.getString("mpi_comm"));
    mesh.metadata().set<size_t>("nb_nodes_including_halo[0]", nodes.size());
    nodes.metadata().set<size_t>("NbRealPts", size_t(nnodes));
    nodes.metadata().set<size_t>("NbVirtualPts", size_t(0));
    nodes.global_index().metadata().set("human_readable", true);
    nodes.global_index().metadata().set("min", 1);
    nodes.global_index().metadata().set("max", nb_nodes_ + grid.ny() + 2);
    mesh.cells().global_index().metadata().set("human_readable", true);
    mesh.cells().global_index().metadata().set("min", 1);
    mesh.cells().global_index().metadata().set("max", nb_points_);

    //generateGlobalElementNumbering(mesh);
}  // generate_mesh

namespace {
static MeshGeneratorBuilder<HealpixMeshGenerator> __HealpixMeshGenerator(HealpixMeshGenerator::static_type());
}

}  // namespace meshgenerator
}  // namespace atlas
