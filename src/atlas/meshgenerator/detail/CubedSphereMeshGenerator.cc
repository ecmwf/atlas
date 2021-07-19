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

  // Set grid shaped arrays of cell local indices.
  array::ArrayT<int> cellLocalIdxData(nTiles, N, N);
  auto cellLocalIdx = array::make_view<int, 3>(cellLocalIdxData);

  // Loop over all cells and set local index
  int iCell = 0;
  for (idx_t t = 0; t < nTiles; ++t) {
    for (idx_t j = 0; j < N; ++j) {
      for (idx_t i = 0; i < N; ++i) {
        cellLocalIdx(t, j, i) = iCell++;
      }
    }
  }

  // Loop over grid points and set cell global index and xy
  idx_t iGrid = 0;
  auto tji = csgrid.tij().begin();
  for (const auto& xy : csgrid.xy()) {

      // Get t, j, i.
      const auto t = (*tji).t(), j = (*tji).j(), i = (*tji).i();
      ++tji;

      // Get local index.
      const auto iLocal = cellLocalIdx(t, j, i);

      // Set cell-centroid xy.
      cellXY[static_cast<size_t>(iLocal)] = xy;

      // Set remote id to -1 (all cells owned for now)
      cellRemoteIdx(iLocal) = -1;

      // Set partition using grid global index.
      cellPart(iLocal) = distribution.partition(iGrid);

      // Set cell global index to grid global index.
      cellGlobalIdx(iLocal) = iGrid++;

  }

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

  // Set grid shaped array of node local indices.
  array::ArrayT<int> nodeLocalIdxData(nTiles, N + 1, N + 1);
  auto nodeLocalIdx = array::make_view<int, 3>(nodeLocalIdxData);

  // Loop over all nodes and set local index
  int iNode = 0;
  for (idx_t t = 0; t < nTiles; ++t) {
    for (idx_t j = 0; j < N + 1; ++j) {
      for (idx_t i = 0; i < N + 1; ++i) {
        nodeLocalIdx(t, j, i) = iNode++;
      }
    }
  }

  // Initialise node connectivity.
  auto& nodeConnectivity = mesh.cells().node_connectivity();

  // Loop over cells and add nodes.
  idx_t iNodeGlobal = 0;
  for (idx_t t = 0; t < nTiles; ++t) {

    // Calc coordinate xy Jacobian wrt ij (this belongs in tile/grid class).
    const auto i0_j0 = static_cast<size_t>(cellLocalIdx(t, 0, 0));
    const auto i1_j0 = static_cast<size_t>(cellLocalIdx(t, 0, 1));
    const auto i0_j1 = static_cast<size_t>(cellLocalIdx(t, 1, 0));

    const auto dx_di = cellXY[i1_j0].x() - cellXY[i0_j0].x();
    const auto dx_dj = cellXY[i0_j1].x() - cellXY[i0_j0].x();
    const auto dy_di = cellXY[i1_j0].y() - cellXY[i0_j0].y();
    const auto dy_dj = cellXY[i0_j1].y() - cellXY[i0_j0].y();

    for (idx_t j = 0; j < N; ++j) {
      for (idx_t i = 0; i < N; ++i) {

        // Node layout relative to cell.
        //
        //  ^  N3 - N2
        //  ^  |  C  |
        //  j  N0 - N1
        //      i > >
        //
        // C: Cell
        // N0, N1, N2, N3: Nodes

        // Get cell.
        const auto iCell = cellLocalIdx(t, j, i);

        const auto iNodes = std::array<idx_t, 4> {
          nodeLocalIdx(t, j, i),
          nodeLocalIdx(t, j, i + 1),
          nodeLocalIdx(t, j + 1, i + 1),
          nodeLocalIdx(t, j + 1, i)
        };

        auto addNode = [&](const idx_t iNode, const idx_t iCell,
          const double di, const double dj) {

          // Get cell centre.
          const auto x0 = cellXY[static_cast<size_t>(iCell)].x();
          const auto y0 = cellXY[static_cast<size_t>(iCell)].y();

          // Set node position.
          nodeXY(iNode, XX) = x0 + di * dx_di + dj * dx_dj;
          nodeXY(iNode, YY) = y0 + di * dy_di + dj * dy_dj;

          // Add some lon-lat cleverness here.

          // Set remote id to -1 (all node owned for now).
          nodeRemoteIdx(iNode) = -1;

          // Set partition to same as cell (for now).
          nodePart(iNode) = cellPart(iCell);

          // Set unique node global index.
          nodeGlobalIdx(iNode) = iNodeGlobal++;

          return;
        };


        // Every cell sets top right node.
        addNode(iNodes[2], iCell, 0.5, 0.5);

        // i == 0 cells set top left node.
        if (i == 0) addNode(iNodes[3], iCell, -0.5, 0.5);

        // j == 0 cells set bottom right node.
        if (j == 0) addNode(iNodes[1], iCell, 0.5, -0.5);

        // i == 0 and j == 0 cell sets bottom left node.
        if (i == 0 && j == 0) addNode(iNodes[0], iCell, -0.5, -0.5);

        // Set node connectivity.
        nodeConnectivity.set(iCell, iNodes.data());

      }
    }
  }

  return;
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
