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
    configure_defaults();

    size_t nb_parts;
    if (p.get("nb_parts", nb_parts)) {
        options.set("nb_parts", nb_parts);
    }

    size_t part;
    if (p.get("part", part)) {
        options.set("part", part);
    }

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
    // This option sets number of parts the mesh will be split in
    options.set("nb_parts", mpi::size());

    // This option sets the part that will be generated
    options.set("part", mpi::rank());

    // This options sets the default partitioner
    std::string partitioner;
    if (grid::Partitioner::exists("equal_regions") && mpi::size() > 1) {
        partitioner = "equal_regions";
    }
    else {
        partitioner = "serial";
    }
    options.set<std::string>("partitioner", partitioner);
}

namespace {
int idx_xy_to_x(const int xidx, const int yidx, const int ns) {
    ATLAS_ASSERT(yidx < 4 * ns + 1 && yidx >= 0);
    ATLAS_ASSERT(xidx >= 0);

    auto ghostIdx = [ns](int latid) { return 12 * ns * ns + 16 + latid; };

    if (yidx == 0) {
        ATLAS_ASSERT(xidx < 9 && xidx >= 0);
        return (xidx != 8 ? xidx : ghostIdx(yidx));
    }
    else if (yidx < ns) {
        ATLAS_ASSERT(xidx < 4 * yidx + 1 && xidx >= 0);
        return (xidx != 4 * yidx ? 2 * yidx * (yidx - 1) + 8 + xidx : ghostIdx(yidx));
    }
    else if (yidx <= 2 * ns) {
        ATLAS_ASSERT(xidx < 4 * ns + 1 && xidx >= 0);
        return (xidx != 4 * ns ? 2 * ns * (ns - 1) + 4 * ns * (yidx - ns) + 8 + xidx : ghostIdx(yidx));
    }
    else if (yidx <= 3 * ns) {
        ATLAS_ASSERT(xidx < 4 * ns + 1 && xidx >= 0);
        return (xidx != 4 * ns ? 2 * ns * (3 * ns + 1) + 4 * ns * (yidx - 2 * ns - 1) + 8 + xidx : ghostIdx(yidx));
    }
    else if (yidx == 3 * ns + 1 && ns > 1) {
        ATLAS_ASSERT(xidx < 4 * (ns - 1) + 1 && xidx >= 0);
        return (xidx != 4 * (ns - 1) ? 2 * ns * (5 * ns + 1) + 4 * ns * (yidx - 3 * ns - 1) + 8 + xidx
                                     : ghostIdx(yidx));
    }
    else if (yidx < 4 * ns) {
        ATLAS_ASSERT(xidx < 4 * (ns - (yidx - 3 * ns)) + 1 && xidx >= 0);
        return (xidx != 4 * (ns - (yidx - 3 * ns)) ? 2 * ns * (5 * ns + 1) + 4 * ns * (yidx - 3 * ns - 1) -
                                                         2 * (yidx - 3 * ns) * (yidx - 3 * ns - 1) + 8 + xidx
                                                   : ghostIdx(yidx));
    }
    else {
        ATLAS_ASSERT(xidx < 9 && xidx >= 0);
        return (xidx != 8 ? 12 * ns * ns + 8 + xidx : ghostIdx(yidx));
    }
}

int up_idx(const int xidx, const int yidx, const int ns) {
    ATLAS_ASSERT(yidx <= 4 * ns && yidx >= 0);

    auto ghostIdx = [ns](int latid) { return 12 * ns * ns + 16 + latid; };

    int ret;

    // global idx
    if (yidx == 0) {
        ATLAS_ASSERT(xidx < 8);
        ret = (xidx != 7 ? xidx + 1 : ghostIdx(0));
    }
    else if (yidx == 1) {
        ATLAS_ASSERT(xidx < 4);
        ret = 2 * xidx;
    }
    else if (yidx < ns) {
        ATLAS_ASSERT(xidx < 4 * yidx);
        if (xidx != 4 * yidx - 1) {
            ret = 2 * (yidx - 2) * (yidx - 1) + 8 + xidx - std::floor(xidx / (double)yidx);
        }
        else {
            ret = ghostIdx(yidx - 1);
        }
    }
    else if (yidx == ns && ns < 3) {
        ATLAS_ASSERT(xidx < 4 * ns);
        if (xidx != 4 * ns - 1) {
            ret = 2 * ns * (ns - 1) + 8 - 4 * (ns - 1) + (xidx + 1) / 2;
        }
        else {
            ret = ghostIdx(yidx - 1);
        }
    }
    else if (yidx == ns) {
        ATLAS_ASSERT(xidx < 4 * ns);
        if (xidx != 4 * ns - 1) {
            ret = 2 * (ns - 2) * (ns - 1) + 8 + xidx - std::floor(xidx / (double)yidx);
        }
        else {
            ret = ghostIdx(yidx - 1);
        }
    }
    else if (yidx <= 3 * ns) {
        ATLAS_ASSERT(xidx < 4 * ns);
        int stg = (yidx - ns) % 2;
        if (xidx != 4 * ns - 1 || (xidx == 4 * ns - 1 && stg)) {
            ret = 2 * ns * (ns - 1) + 8 + 4 * ns * (yidx - ns - 1) + xidx + 1 - stg;
        }
        else {
            ret = ghostIdx(yidx - 1);
        }
    }
    else if (yidx < 4 * ns - 1) {
        int yidxl = 4 * ns - yidx;
        ATLAS_ASSERT(xidx < 4 * yidxl);
        ret = 12 * ns * ns + 9 - 2 * (yidxl + 2) * (yidxl + 1) + xidx + std::floor(xidx / (double)yidxl);
    }
    else if (yidx == 4 * ns - 1) {
        ATLAS_ASSERT(xidx < 4);
        ret = 12 * ns * ns + 5 - (ns == 1 ? 4 : 8) + 2 * xidx;
    }
    else {
        ATLAS_ASSERT(xidx < 8);
        if (ns == 1) {
            if (xidx != 7) {
                ret = 12 * ns * ns + 4 + (xidx % 2 ? -4 + (xidx + 1) / 2 : xidx / 2);
            }
            else {
                ret = ghostIdx(4 * ns - 2);
            }
        }
        else {
            ret = 12 * ns * ns + 8 + (xidx % 2 ? xidx - 12 : xidx - 4 - xidx / 2);
        }
    }
    return ret;
}

int down_idx(const int xidx, const int yidx, const int ns) {
    ATLAS_ASSERT(yidx <= 4 * ns);

    auto ghostIdx = [ns](int latid) { return 12 * ns * ns + 16 + latid; };

    int ret;

    // global idx
    if (yidx == 0) {
        if (xidx < 8) {
            ATLAS_ASSERT(xidx < 8);
        }
        if (ns == 1) {
            if (xidx != 7) {
                ret = 8 + (xidx % 2 ? 4 + (xidx + 1) / 2 : xidx / 2);
            }
            else {
                ret = ghostIdx(2);
            }
        }
        else {
            ret = 8 + ((xidx + 1) % 2 ? xidx / 2 : 4 + xidx);
        }
    }
    else if (yidx < ns) {
        ATLAS_ASSERT(xidx < 4 * yidx);
        ret = 2 * yidx * (yidx + 1) + 9 + xidx + std::floor(xidx / (double)yidx);
    }
    else if (yidx == ns && ns == 1) {
        ATLAS_ASSERT(xidx < 4);
        ret = (xidx != 3 ? 13 + xidx : ghostIdx(2));
    }
    else if (yidx == 2 * ns && ns == 1) {
        ATLAS_ASSERT(xidx < 4);
        ret = 16 + xidx;
    }
    else if (yidx < 3 * ns && ns > 1) {
        ATLAS_ASSERT(xidx < 4 * ns);
        int stg = (yidx - ns) % 2;
        if (xidx != 4 * ns - 1 || (xidx == 4 * ns - 1 && stg)) {
            ret = 2 * ns * (ns - 1) + 8 + 4 * ns * (yidx - ns + 1) + xidx + (yidx != 3 * ns ? 1 - stg : 0);
        }
        else {
            ret = ghostIdx(yidx + 1);
        }
    }
    else if (yidx == 4 * ns - 2) {
        ATLAS_ASSERT(xidx < 8);
        ret = (xidx != 7 ? 12 * ns * ns + 4 + (xidx + 1) / 2 : ghostIdx(4 * ns - 1));
    }
    else if (yidx == 4 * ns - 1) {
        ATLAS_ASSERT(xidx < 4);
        ret = 12 * ns * ns + 8 + 2 * xidx;
    }
    else if (yidx < 4 * ns - 1) {
        int yidxl = yidx - 3 * ns;
        ATLAS_ASSERT(xidx < 4 * (ns - yidxl));
        if (xidx != 4 * (ns - yidxl) - 1) {
            ret = 2 * ns * (5 * ns + 1) + 8 + 4 * ns * yidxl - 2 * (yidxl + 1) * yidxl + xidx -
                  std::floor(xidx / (double)(ns - yidxl));
        }
        else {
            ret = ghostIdx(yidx + 1);
        }
    }
    else if (yidx == 4 * ns) {
        ATLAS_ASSERT(xidx < 8);
        ret = (xidx != 7 ? 12 * ns * ns + 8 + xidx + 1 : ghostIdx(yidx));
    }
    else {
        throw_AssertionFailed("Invalid value of yidx", Here());
    }
    return ret;
}

int right_idx(const int xidx, const int yidx, const int ns) {
    ATLAS_ASSERT(yidx <= 4 * ns);

    auto ghostIdx = [ns](int latid) { return 12 * ns * ns + 16 + latid; };
    int ret       = -1;

    if (yidx == 0) {
        if (xidx < 8) {
            ATLAS_ASSERT(xidx < 8);
        }
        if (ns == 1) {
            ret = (xidx != 7 ? (xidx % 2 ? 8 + (xidx + 1) / 2 : 13 + xidx / 2) : ghostIdx(1));
        }
        else {
            ret = (xidx < 7 ? (xidx % 2 ? 8 + (xidx + 1) / 2 : 13 + xidx) : ghostIdx(1));
        }
    }
    else if (yidx == 1) {
        ATLAS_ASSERT(xidx < 4);
        ret = (xidx < 7 ? 1 + 2 * xidx : ghostIdx(0));
    }
    else if (yidx < ns) {
        ATLAS_ASSERT(xidx < 4 * yidx);
        ret = (xidx != 4 * yidx - 1 ? 2 * yidx * (yidx - 1) + 9 + xidx : ghostIdx(yidx));
    }
    else if (yidx == 3 && ns == 1) {
        ATLAS_ASSERT(xidx < 4);
        ret = 21 + 2 * xidx;
    }
    else if (yidx <= 3 * ns) {
        ATLAS_ASSERT(xidx < 4 * ns + 1);
        ret = (xidx != 4 * ns - 1 ? 2 * ns * (ns - 1) + 4 * ns * (yidx - ns) + 9 + xidx : ghostIdx(yidx));
    }
    else if (yidx < 4 * ns - 1 && ns > 1) {
        int yidxl = yidx - 3 * ns;
        ATLAS_ASSERT(xidx < 4 * (ns - yidxl));
        if (xidx != 4 * (ns - yidxl) - 1) {
            ret = 2 * ns * (5 * ns + 1) + 4 * ns * (yidx - 3 * ns - 1) - 2 * (yidx - 3 * ns) * (yidx - 3 * ns - 1) + 9 +
                  xidx;
        }
        else {
            ret = ghostIdx(yidx);
        }
    }
    else if (yidx == 4 * ns - 1) {
        ATLAS_ASSERT(xidx < 4);
        ret = 12 * ns * ns + 9 + 2 * xidx;
    }
    else if (yidx == 4 * ns) {
        ATLAS_ASSERT(xidx < 8);
        if (xidx != 7) {
            ret = (xidx % 2 ? 12 * ns * ns + 4 + (xidx + 1) / 2 : 12 * ns * ns + 4 - (ns == 1 ? 3 : 7) + xidx);
        }
        else {
            ret = ghostIdx(yidx - 1);
        }
    }
    return ret;
}
}  // namespace

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

    grid::Partitioner partitioner(partitioner_type, nb_parts);
    grid::Distribution distribution(partitioner.partition(grid));

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
    //      mesh.cells().add( new mesh::temporary::Quadrilateral(), nquads  );
    //    further define cells with
    //      array::ArrayView<gidx_t,1> cells_glb_idx( mesh.cells().global_index()
    //      );
    //      array::ArrayView<int,1>    cells_part(    mesh.cells().partition() );
    // - define connectivity with
    //      mesh::HybridElements::Connectivity& node_connectivity =
    //      mesh.cells().node_connectivity();
    //      node_connectivity.set( jcell, quad_nodes );
    //    where quad_nodes is a 4-element integer array containing the LOCAL
    //    indices of the nodes
    //
    // The rule do determine if a cell belongs to a proc is the following with some
    // exceptions near poles: cell belongs to a proc if the left corner of the cell belongs to that proc.


    ATLAS_TRACE();

    ATLAS_ASSERT(HealpixGrid(grid));

    const int mypart    = options.get<size_t>("part");
    const int nparts    = options.get<size_t>("nb_parts");
    const int ny        = grid.ny() + 2;
    const int ns        = (ny - 1) / 4;
    const int nvertices = 12 * ns * ns + 16;

    int inode;
    auto latPoints = [ny, &grid](int latid) { return (latid == 0 ? 8 : (latid == ny - 1 ? 8 : grid.nx()[latid - 1])); };

    int ii, ix, iy, ii_ghost, ii_glb;
    int iy_min, iy_max;   // a belt (iy_min:iy_max) surrounding the nodes on this processor
    int nnodes_nonghost;  // non-ghost node: belongs to this part

    // vector of local indices: necessary for remote indices of ghost nodes
    std::vector<int> local_idx(nvertices, -1);
    std::vector<int> current_idx(nparts, 0);  // index counter for each proc

    // loop over all points to determine local indices and surrounding rectangle
    ii_glb          = 0;
    iy_min          = ny + 1;
    iy_max          = 0;
    nnodes_nonghost = 0;
    for (iy = 0; iy < ny; iy++) {
        int nx = latPoints(iy);
        for (ix = 0; ix < nx; ix++) {
            int proc_id = (iy == 0 ? 0 : (iy == ny - 1 ? mpi::comm().size() - 1 : distribution.partition(ii_glb - 8)));
            local_idx[ii_glb] = current_idx[proc_id]++;
            if (proc_id == mypart) {
                ++nnodes_nonghost;
                iy_min = std::min(iy_min, iy);
                iy_max = std::max(iy_max, iy);
            }
            ++ii_glb;  // global index
        }
    }

#if DEBUG_OUTPUT_DETAIL
    inode = 0;
    Log::info() << "local_idx : " << std::endl;
    for (size_t ilat = 0; ilat < ny; ilat++) {
        for (size_t ilon = 0; ilon < latPoints(ilat); ilon++) {
            Log::info() << std::setw(4) << local_idx[inode];
            inode++;
        }
        Log::info() << std::endl;
    }
    inode = 0;
    Log::info() << "global_idx : " << std::endl;
    for (size_t ilat = 0; ilat < ny; ilat++) {
        for (size_t ilon = 0; ilon < latPoints(ilat); ilon++) {
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
        nnodes_SB += latPoints(iy) + 1;
    }

#if DEBUG_OUTPUT
    std::cout << "[" << mypart << "] : nnodes_SB = " << nnodes_SB << "\n";
#endif

    // partitions and local indices in SB
    std::vector<int> parts_SB(nnodes_SB, -1);
    std::vector<int> local_idx_SB(nnodes_SB, -1);
    std::vector<bool> is_ghost_SB(nnodes_SB, true);

    // global starting node index for the partition
    int parts_sidx = idx_xy_to_x(0, iy_min, ns);

    auto compute_part = [&](int ix, int iy, gidx_t ii_glb) -> int {
        if (ii_glb < 8) {
            // HACK! expects equal_regions partitioner. Better would be partition of attached element of which this node would be the North-West point.
            return 0;
        }
        if (ii_glb > nvertices - 9) {
            // HACK! expects equal_regions partitioner. Better would be partition of attached element of which this node would be the South-West point.
            // Also, we should not have mpi here.
            return mpi::comm().size() - 1;
        }
        return distribution.partition(idx_xy_to_x(ix, iy, ns) - 8);
    };

    ii       = 0;  // index inside SB
    ii_ghost = nnodes_SB - (iy_max - iy_min + 1);
    for (iy = iy_min; iy <= iy_max; iy++) {
        int nx    = latPoints(iy) + 1;
        ii_glb    = ii + parts_sidx;
        int part0 = compute_part(0, iy, ii_glb);
        for (ix = 0; ix < nx; ix++) {
            if (ix != nx - 1) {
                ii_glb           = ii + parts_sidx;
                parts_SB[ii]     = compute_part(ix, iy, ii_glb);
                local_idx_SB[ii] = ii;
                is_ghost_SB[ii]  = !((parts_SB[ii] == mypart));
                ++ii;
            }
            else {
                parts_SB[ii_ghost]     = part0;
                local_idx_SB[ii_ghost] = ii_ghost;
                is_ghost_SB[ii_ghost]  = true;
                ++ii_ghost;
            }
        }
    }

#if DEBUG_OUTPUT_DETAIL
    std::cout << "[" << mypart << "] : "
              << "parts_SB = ";
    for (ii = 0; ii < nnodes_SB; ii++) {
        std::cout << parts_SB[ii] << ",";
    }
    std::cout << std::endl;
    std::cout << "[" << mypart << "] : "
              << "local_idx_SB = ";
    for (ii = 0; ii < nnodes_SB; ii++) {
        std::cout << local_idx_SB[ii] << ",";
    }
    std::cout << std::endl;
    std::cout << "[" << mypart << "] : "
              << "is_ghost_SB = ";
    for (ii = 0; ii < nnodes_SB; ii++) {
        std::cout << is_ghost_SB[ii] << ",";
    }
    std::cout << std::endl;
#endif

    // vectors marking nodes that are necessary for this proc's cells
    std::vector<bool> is_node_SB(nnodes_SB, false);

    // determine number of cells and number of nodes
    int nnodes = 0;
    int ncells = 0;
    ii         = 0;
    int iil;
    for (iy = iy_min; iy <= iy_max; iy++) {
        int nx = latPoints(iy);
        for (ix = 0; ix < nx; ix++) {
            int is_cell = (iy == 0 ? ix % 2 : 1) * (iy == ny - 1 ? ix % 2 : 1);

            if (!is_ghost_SB[ii] && is_cell) {
                // mark this node as being used
                if (!is_node_SB[ii]) {
                    ++nnodes;
                    is_node_SB[ii] = true;
                }

                ++ncells;

                const int glb2loc_ghost_offset = -nnodes_SB + iy_max + 12 * ns * ns + 17;
                // mark upper corner
                iil = up_idx(ix, iy, ns);
                iil -= (iil < 12 * ns * ns + 16 ? parts_sidx : glb2loc_ghost_offset);
                if (!is_node_SB[iil]) {
                    ++nnodes;
                    is_node_SB[iil] = true;
                }
                // mark lower corner
                iil = down_idx(ix, iy, ns);
                iil -= (iil < 12 * ns * ns + 16 ? parts_sidx : glb2loc_ghost_offset);
                if (!is_node_SB[iil]) {
                    ++nnodes;
                    is_node_SB[iil] = true;
                }
                // mark right corner
                iil = right_idx(ix, iy, ns);
                iil -= (iil < 12 * ns * ns + 16 ? parts_sidx : glb2loc_ghost_offset);
                if (!is_node_SB[iil]) {
                    ++nnodes;
                    is_node_SB[iil] = true;
                }
            }
            ++ii;
        }
    }

    // periodic points are always needed, even if they don't belong to a cell
    ii_ghost = nnodes_SB - (iy_max - iy_min + 1);
    for (ii = 0; ii < iy_max - iy_min + 1; ii++) {
        if (!is_node_SB[ii_ghost + ii]) {
            // is_node_SB[ii_ghost + ii] = true;
            // ++nnodes;
        }
    }

#if DEBUG_OUTPUT
    std::cout << "[" << mypart << "] : "
              << "nnodes = " << nnodes << ", ncells = " << ncells << ", parts_sidx = " << parts_sidx << std::endl;
#endif
#if DEBUG_OUTPUT_DETAIL
    std::cout << "[" << mypart << "] : "
              << "is_node_SB = ";
    for (int ii = 0; ii < nnodes_SB; ii++) {
        std::cout << is_node_SB[ii] << ",";
    }
    std::cout << std::endl;
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

    // define cells and associated properties
    mesh.cells().add(new mesh::temporary::Quadrilateral(), ncells);
    int quad_begin                                        = mesh.cells().elements(0).begin();
    auto cells_part                                       = array::make_view<int, 1>(mesh.cells().partition());
    mesh::HybridElements::Connectivity& node_connectivity = mesh.cells().node_connectivity();

    idx_t quad_nodes[4];
    int jcell = quad_begin;

    int inode_nonghost, inode_ghost;

    // loop over nodes and set properties
    inode_nonghost = 0;
    inode_ghost    = nnodes_nonghost;  // ghost nodes start counting after nonghost nodes
    ii             = 0;                // index inside SB
    for (iy = iy_min; iy <= iy_max; iy++) {
        int nx = latPoints(iy) + 1;
        for (ix = 0; ix < nx; ix++) {
            int iil = idx_xy_to_x(ix, iy, ns);
            iil -= (iil < 12 * ns * ns + 16 ? parts_sidx : -nnodes_SB + iy_max + 12 * ns * ns + 17);
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

                glb_idx(inode) = idx_xy_to_x(ix, iy, ns) + 1;

                // grid coordinates
                double _xy[2];
                double xy1[2], xy2[2];
                if (iy == 0) {
                    _xy[0] = 45. * ix;
                    _xy[1] = 90.;
                    Topology::set(flags(inode), Topology::BC);
                }
                else if (iy == ny - 1) {
                    _xy[0] = 45. * ix;
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
                std::cout << "[" << mypart << "] : "
                          << "New node \tinode=" << inode << "; iil= " << iil << "; ix=" << ix << "; iy=" << iy
                          << "; glon=" << lonlat(inode, 0) << "; glat=" << lonlat(inode, 1)
                          << "; glb_idx=" << glb_idx(inode) << "; loc_idx=" << local_idx_SB[iil] << std::endl;
#endif
            }
            ii += (ix != nx - 1 ? 1 : 0);
        }
    }

    ii = 0;  // index inside SB (surrounding belt)
    for (iy = iy_min; iy <= iy_max; iy++) {
        int nx = latPoints(iy) + 1;
        for (ix = 0; ix < nx; ix++) {
            int is_cell = (iy == 0 ? ix % 2 : 1) * (iy == ny - 1 ? ix % 2 : 1);
            int iil     = idx_xy_to_x(ix, iy, ns);
            iil -= (iil < 12 * ns * ns + 16 ? parts_sidx : -nnodes_SB + iy_max + 12 * ns * ns + 17);
            if (!is_ghost_SB[iil] && is_cell) {
                // define cell corners (local indices)
                quad_nodes[0] = local_idx_SB[iil];

                quad_nodes[1] = down_idx(ix, iy, ns);  // point to the right
                quad_nodes[1] -=
                    (quad_nodes[1] < 12 * ns * ns + 16 ? parts_sidx : -nnodes_SB + iy_max + 12 * ns * ns + 17);
                quad_nodes[1] = local_idx_SB[quad_nodes[1]];

                quad_nodes[2] = right_idx(ix, iy, ns);  // point above right
                quad_nodes[2] -=
                    (quad_nodes[2] < 12 * ns * ns + 16 ? parts_sidx : -nnodes_SB + iy_max + 12 * ns * ns + 17);
                quad_nodes[2] = local_idx_SB[quad_nodes[2]];

                quad_nodes[3] = up_idx(ix, iy, ns);  // point above
                quad_nodes[3] -=
                    (quad_nodes[3] < 12 * ns * ns + 16 ? parts_sidx : -nnodes_SB + iy_max + 12 * ns * ns + 17);
                quad_nodes[3] = local_idx_SB[quad_nodes[3]];

                node_connectivity.set(jcell, quad_nodes);
                cells_part(jcell) = mypart;
#if DEBUG_OUTPUT_DETAIL
                std::cout << "[" << mypart << "] : "
                          << "New quad " << jcell << ": " << glb_idx(quad_nodes[0]) << "," << glb_idx(quad_nodes[1])
                          << "," << glb_idx(quad_nodes[2]) << "," << glb_idx(quad_nodes[3]) << std::endl;
#endif
                ++jcell;
            }
            ii += (ix != nx - 1 ? 1 : 0);
        }
    }

#if DEBUG_OUTPUT_DETAIL
    // list nodes
    for (inode = 0; inode < nnodes; inode++) {
        std::cout << "[" << mypart << "] : "
                  << " node " << inode << ": ghost = " << ghost(inode) << ", glb_idx = " << glb_idx(inode)
                  << ", part = " << part(inode) << ", lon = " << lonlat(inode, 0) << ", lat = " << lonlat(inode, 1)
                  << ", remote_idx = " << remote_idx(inode) << std::endl;
    }

    int* cell_nodes;
    for (jcell = 0; jcell < ncells; jcell++) {
        std::cout << "[" << mypart << "] : "
                  << " cell " << jcell << ": " << glb_idx(node_connectivity(jcell, 0)) << ","
                  << glb_idx(node_connectivity(jcell, 1)) << "," << glb_idx(node_connectivity(jcell, 2)) << ","
                  << glb_idx(node_connectivity(jcell, 3)) << std::endl;
    }
#endif

    mesh.metadata().set<size_t>("nb_nodes_including_halo[0]", nodes.size());
    nodes.metadata().set<size_t>("NbRealPts", size_t(nnodes));
    nodes.metadata().set<size_t>("NbVirtualPts", size_t(0));
    nodes.global_index().metadata().set("human_readable", true);
    nodes.global_index().metadata().set("min", 1);
    nodes.global_index().metadata().set("max", nvertices + grid.ny() + 2);


    generateGlobalElementNumbering(mesh);

}  // generate_mesh

namespace {
static MeshGeneratorBuilder<HealpixMeshGenerator> __HealpixMeshGenerator(HealpixMeshGenerator::static_type());
}

}  // namespace meshgenerator
}  // namespace atlas
