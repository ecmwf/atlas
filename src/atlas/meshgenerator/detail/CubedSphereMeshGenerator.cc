/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <vector>

#include "eckit/utils/Hash.h"

#include "atlas/field/Field.h"
#include "atlas/grid/CubedSphereGrid.h"
#include "atlas/grid/Distribution.h"
#include "atlas/grid/Iterator.h"
#include "atlas/grid/Partitioner.h"
#include "atlas/library/config.h"
#include "atlas/mesh/Connectivity.h"
#include "atlas/mesh/ElementType.h"
#include "atlas/mesh/Elements.h"
#include "atlas/mesh/HybridElements.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/meshgenerator/detail/CubedSphereMeshGenerator.h"
#include "atlas/meshgenerator/detail/MeshGeneratorFactory.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"
#include "atlas/util/Constants.h"
#include "atlas/util/CoordinateEnums.h"

#define DEBUG_OUTPUT 0
#define DEBUG_OUTPUT_DETAIL 0

using Topology = atlas::mesh::Nodes::Topology;

namespace atlas {
namespace meshgenerator {

// -------------------------------------------------------------------------------------------------

CubedSphereMeshGenerator::CubedSphereMeshGenerator(const eckit::Parametrisation& p) {}

// -------------------------------------------------------------------------------------------------

void CubedSphereMeshGenerator::configure_defaults() {}

// -------------------------------------------------------------------------------------------------

void CubedSphereMeshGenerator::generate(const Grid& grid, Mesh& mesh) const {

  // Check for proper grid and need for mesh
  ATLAS_ASSERT(!mesh.generated());
  if (!CubedSphereGrid(grid)) {
    throw_Exception("CubedSphereMeshGenerator can only work "
      "with a cubedsphere grid", Here());
  }

  // Check for proper stagger
  const auto gridType = grid->type();
  const auto gridStagger = gridType.substr(gridType.rfind("-") - 1, 1);

  if (gridStagger != "C") {
    throw_Exception("CubedSphereMeshGenerator can only work with a"
      "cell-centroid grid. Try FV3CubedSphereMeshGenerator instead.");
  }

  // Partitioner
  const auto partitioner = grid::Partitioner("checkerboard", 1);
  const auto distribution =
    grid::Distribution(partitioner.partition(grid));

  generate(grid, distribution, mesh);
}

// -------------------------------------------------------------------------------------------------

void CubedSphereMeshGenerator::generate(const Grid& grid, const grid::Distribution& distribution, Mesh& mesh) const {
    const auto csgrid = CubedSphereGrid(grid);

    // Get dimensions of grid
    const int N      = csgrid.N();
    const int nTiles = csgrid.GetNTiles();

    ATLAS_TRACE("CubedSphereMeshGenerator::generate");
    Log::debug() << "Number of faces per tile edge = " << std::to_string(N) << std::endl;

    // Number of nodes and cells.
    const int nNodesUnique  = nTiles * N * N + 2;
    const int nNodesAll     = nTiles * (N + 1) * (N + 1);
    const int nCells        = nTiles * N * N;

    // Construct mesh cells.
    auto& cells = mesh.cells();
    cells.add(new mesh::temporary::Quadrilateral(), nCells);

    auto cellGlobalIdx = array::make_indexview<idx_t, 1>(cells.global_index());
    auto cellRemoteIdx = array::make_indexview<idx_t, 1>(cells.remote_index());
    auto cellPart      = array::make_view<int, 1>(cells.partition());
    auto cellXY        = std::vector<PointXY>(static_cast<size_t>(nCells));

    // Construct mesh nodes.
    auto& nodes = mesh.nodes();
    nodes.resize(nNodesAll);

    auto nodeGlobalIdx = array::make_indexview<idx_t, 1>(nodes.global_index());
    auto nodeRemoteIdx = array::make_indexview<idx_t, 1>(nodes.remote_index());
    auto nodePart      = array::make_view<int, 1>(nodes.partition());
    auto nodeXY        = array::make_view<double, 2>(nodes.xy());
    auto nodeLonLat    = array::make_view<double, 2>(nodes.lonlat());
    auto nodeGhost     = array::make_view<int, 1>(nodes.ghost());
    auto nodeFlags     = array::make_view<int, 1>(nodes.flags());

    // Set grid shaped arrays of cell and node local index
    array::ArrayT<int> cellLocalIdxData(nTiles, N, N);
    array::ArrayT<int> nodeLocalIdxData(nTiles, N + 1, N + 1);
    auto cellLocalIdx = array::make_view<int, 3>(cellLocalIdxData);
    auto nodeLocalIdx = array::make_view<int, 3>(nodeLocalIdxData);



    // Loop over all cells and set local index
    int iCell = 0;
    for (idx_t t = 0; t < nTiles; ++t) {
      for (idx_t j = 0; j < N; ++j) {
        for (idx_t i = 0; i < N; ++i) {
          cellLocalIdx(t, j, i) = iCell++;
        }
      }
    }

    // Loop over all nodes and set local index
    int iNode = 0;
    for (idx_t t = 0; t < nTiles; ++t) {
      for (idx_t j = 0; j < N + 1; ++j) {
        for (idx_t i = 0; i < N + 1; ++i) {
          nodeLocalIdx(t, j, i) = iNode++;
        }
      }
    }

    // Loop over grid points and set cell global index and xy
    idx_t iGrid = 0;
    auto tji = csgrid.tij().begin();
    for (auto& xy : csgrid.xy()) {

      const int iLocalIdx = cellLocalIdx(*(tji).t(), tij.j(), tij.i());

    }





    for (auto& tij : csgrid.tij()) {

      const int iLocalIdx = cellLocalIdx(tij.t(), tij.j(), tij.i());

      cellGlobalIdx(iLocalIdx) = iGrid++;
      cellXY[static_cast<size_t>(iLocalIdx)] =
        csgrid.xy(tij.i(), tij.j(), tij.t());

    }

    // Loop over cells and add nodes.
    for (idx_t t = 0; t < nTiles; ++t) {

      for (idx_t j = 0; j < N; ++j) {
        for (idx_t i = 0; i < N; ++i) {




        }
      }
    }















    // Loop over entire grid
    // ---------------------

    // Node array will include duplicates of the grid points to make it easier to fill up the
    // neighbours. However the mesh will not contain these points so no Ghost nodes are required.
    // We could include the duplicate points in the array if we make them Ghost nodes but it is not
    // clear if this provides any benefit.

    array::ArrayT<int> NodeArrayT(nTiles, N + 1, N + 1);  // All grid points including duplicates
    auto NodeArray = array::make_view<int, 3>(NodeArrayT);

    // Add owned nodes to node array
    idx_t nOwned = 0;
    auto addOwnedNode = [&](Tij tijOwned) {

        // Set node array
        NodeArray(tijOwned[0], tijOwned[1], tijOwned[2]) = nOwned;

        // Get xy from global xy grid array
        double xy_[3];
        csgrid.xy(tijOwned[1], tijOwned[2], tijOwned[0], xy_);

        xy(nOwned, XX) = xy_[XX];
        xy(nOwned, YY) = xy_[YY];

        // Get lonlat from global lonlat array
        double lonlat_[2];
        csgrid.lonlat(tijOwned[1], tijOwned[2], tijOwned[0], lonlat_);

        lonlat(nOwned, LON) = lonlat_[LON];
        lonlat(nOwned, LAT) = lonlat_[LAT];

        // Is not ghost node
        mesh::Nodes::Topology::reset(flags(nOwned));
        ghost(nOwned) = 0;

        glb_idx(nOwned) = nOwned + 1;
        remote_idx(nOwned) = nOwned;
        part(nOwned) = distribution.partition(nOwned);

        ++nOwned;

        return;
    };

    // Loop over owned (t, i, j)

    for (auto& p : csgrid.tij()) addOwnedNode(Tij{p.t(), p.i(), p.j()});

    // Assert that the correct number of nodes have been set
    ATLAS_ASSERT(nnodes == nOwned, "Insufficient nodes");

    // Vector of ghost global index of each ghost point
    auto ghostGblIdx = std::vector<idx_t>();
    auto ownedGblIdx = std::vector<idx_t>();

    // Add ghost nodes to node array
    // (nGhost started after nOwned)
    auto nGhost = nOwned;
    auto addGhostNode = [&](Tij tijGhost, Tij tijOwned) {

        // Get concrete node id
        auto nOwned = NodeArray(tijOwned[0], tijOwned[1], tijOwned[2]);

        // Add ghost node to NodeArray
        NodeArray(tijGhost[0], tijGhost[1], tijGhost[2]) = nGhost;

        // "Create" ghost xy coordinate.

        // Get Jacobian of coords rtw indices
        auto x0 = xy(NodeArray(tijGhost[0], 0, 0), XX);
        auto y0 = xy(NodeArray(tijGhost[0], 0, 0), YY);

        auto dx_di = xy(NodeArray(tijGhost[0], 1, 0), XX) - x0;
        auto dx_dj = xy(NodeArray(tijGhost[0], 0, 1), XX) - x0;
        auto dy_di = xy(NodeArray(tijGhost[0], 1, 0), YY) - y0;
        auto dy_dj = xy(NodeArray(tijGhost[0], 0, 1), YY) - y0;

        // Set xy coordinates
        xy(nGhost, XX) = x0 + tijGhost[1] * dx_di + tijGhost[2] * dx_dj;
        xy(nGhost, YY) = y0 + tijGhost[1] * dy_di + tijGhost[2] * dy_dj;

        // Same lonlat as concrete points
        lonlat(nGhost, LON) = lonlat(nOwned, LON);
        lonlat(nGhost, LAT) = lonlat(nOwned, LAT);

        // Is ghost node
        mesh::Nodes::Topology::set(flags(nGhost), mesh::Nodes::Topology::GHOST);
        ghost(nGhost) = 1;

        // Partitioning logic to be added in future PR
        glb_idx(nGhost) = nGhost + 1;
        remote_idx(nGhost) = nGhost;

        // not sure (below - for multiple PEs)
        part(nGhost) = distribution.partition(nOwned);

        // Append metadata
        // Global indicies of ghost node and owned node
        ghostGblIdx.push_back(nGhost + 1);
        ownedGblIdx.push_back(nOwned + 1);

        ++nGhost;

    };

    // Loop over ghost (t, i, j)
    for (auto& pPair : ghostToOwnedTij) addGhostNode(pPair.first, pPair.second);


    // Assert that the correct number of nodes have been set when duplicates are added
    ATLAS_ASSERT(nnodes_all == nGhost, "Insufficient nodes");

    for (idx_t it = 0; it < nTiles; it++) {
        for (idx_t ix = 0; ix < N + 1; ix++) {
            for (idx_t iy = 0; iy < N + 1; iy++) {
                ATLAS_ASSERT(NodeArray(it, ix, iy) != -9999, "Node Array Not Set Properly");
            }
        }
    }

    // Cells in mesh
    mesh.cells().add(new mesh::temporary::Quadrilateral(), nTiles * N * N);
    //int quad_begin  = mesh.cells().elements(0).begin();
    auto cells_part = array::make_view<int, 1>(mesh.cells().partition());
    auto cells_gidx = array::make_view<gidx_t, 1>(mesh.cells().global_index());
    auto cells_ridx = array::make_indexview<idx_t, 1>(mesh.cells().remote_index());
    atlas::mesh::HybridElements::Connectivity& node_connectivity = mesh.cells().node_connectivity();

    int icell = 0;
    idx_t quad_nodes[4];

    for (int it = 0; it < nTiles; it++) {
        for (int ix = 0; ix < N; ix++) {
            for (int iy = 0; iy < N; iy++) {
                quad_nodes[0] = NodeArray(it, ix, iy);
                quad_nodes[1] = NodeArray(it, ix + 1, iy);
                quad_nodes[2] = NodeArray(it, ix + 1, iy + 1);
                quad_nodes[3] = NodeArray(it, ix, iy + 1);

                node_connectivity.set(icell, quad_nodes);
                cells_part(icell) = part(quad_nodes[0]);
                cells_gidx(icell) = icell + 1;
                cells_ridx(icell) = icell;

                ++icell;
            }
        }
    }

    // Assertion that correct number of cells are set
    ATLAS_ASSERT(ncells == icell, "Insufficient cells have been set");

    // Parallel
    generateGlobalElementNumbering(mesh);
    nodes.metadata().set("parallel", true);

    // Global indices of ghost nodes.
    nodes.metadata().set("ghost-global-idx", ghostGblIdx);

    // Global indices of owned nodes for each ghost node (same order as above)
    nodes.metadata().set("owned-global-idx", ownedGblIdx);
}

// -------------------------------------------------------------------------------------------------

void CubedSphereMeshGenerator::hash(eckit::Hash& h) const {
    h.add("CubedSphereMeshGenerator");
    options.hash(h);
}

// -------------------------------------------------------------------------------------------------

namespace {
static MeshGeneratorBuilder<CubedSphereMeshGenerator> CubedSphereMeshGenerator(
        CubedSphereMeshGenerator::static_type());
}

// -------------------------------------------------------------------------------------------------

}  // namespace meshgenerator
}  // namespace atlas
