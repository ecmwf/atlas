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
#include "atlas/grid/Tiles.h"
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
#include "atlas/projection/detail/CubedSphereProjectionBase.h"
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
  const auto gridType = grid.name();

  std::cout << gridType << std::endl;
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
  const auto csGrid = CubedSphereGrid(grid);

  std::cout << "set dimensions" << std::endl;

  // Get dimensions of grid
  const int N      = csGrid.N();
  const int nTiles = csGrid.GetNTiles();

  // Set projection and tiles.

  std::cout << "get tiles" << std::endl;
  auto gridTiles = CubedSphereTiles(util::Config("type", "cubedsphere_lfric"));

  ATLAS_TRACE("CubedSphereMeshGenerator::generate");
  Log::debug() << "Number of faces per tile edge = " << std::to_string(N) << std::endl;

  // Number of nodes and cells.
  const int nNodesUnique  = nTiles * N * N + 2;
  const int nNodesAll     = nTiles * (N + 1) * (N + 1);
  const int nCells        = nTiles * N * N;

  std::cout << "construct cells" << std::endl;

  // Construct mesh cells.
  auto& cells = mesh.cells();
  cells.add(new mesh::temporary::Quadrilateral(), nCells);

  std::cout << "added quads" << std::endl;

  auto cellRemoteIdx = array::make_indexview<idx_t, 1>(cells.remote_index());
  auto cellGlobalIdx = array::make_view<gidx_t, 1>(cells.global_index());
  auto cellPart      = array::make_view<int, 1>(cells.partition());

  std::cout << "pointed views" << std::endl;

  auto cellXY        = std::vector<PointXY>(static_cast<size_t>(nCells));

   std::cout << "made vectors" << std::endl;

  // Set grid shaped arrays of cell local indices.
  array::ArrayT<int> cellLocalIdxData(nTiles, N, N);
  auto cellLocalIdx = array::make_view<int, 3>(cellLocalIdxData);


  std::cout << "populate cells" << std::endl;

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
  gidx_t iGrid = 1;
  auto tji = csGrid.tij().begin();
  for (const auto& xy : csGrid.xy()) {

      // Get t, j, i.
      const auto t = (*tji).t(), j = (*tji).j(), i = (*tji).i();
      ++tji;

      // Get local index.
      const auto iCellLocal = cellLocalIdx(t, j, i);

      // Set cell-centroid xy.
      cellXY[static_cast<size_t>(iCellLocal)] = xy;

      // Set remote id to iLocal (all cells owned for now)
      cellRemoteIdx(iCellLocal) = iCellLocal;

      // Set partition using grid global index.
      cellPart(iCellLocal) = 0 ;

      // Set cell global index to grid global index.
      cellGlobalIdx(iCellLocal) = iGrid++;

  }

  std::cout << "construct nodes" << std::endl;
  // Construct mesh nodes.
  auto& nodes = mesh.nodes();
  nodes.resize(nNodesAll);

  auto nodeRemoteIdx = array::make_indexview<idx_t, 1>(nodes.remote_index());
  auto nodeGlobalIdx = array::make_view<gidx_t, 1>(nodes.global_index());
  auto nodePart      = array::make_view<int, 1>(nodes.partition());
  auto nodeXY        = array::make_view<double, 2>(nodes.xy());
  auto nodeLonLat    = array::make_view<double, 2>(nodes.lonlat());
  auto nodeGhost     = array::make_view<int, 1>(nodes.ghost());
  auto nodeFlags     = array::make_view<int, 1>(nodes.flags());

  // Set grid shaped array of node local indices.
  array::ArrayT<int> nodeLocalIdxData(nTiles, N + 1, N + 1);
  auto nodeLocalIdx = array::make_view<int, 3>(nodeLocalIdxData);


  // Initialise node connectivity.
  auto& nodeConnectivity = mesh.cells().node_connectivity();

  std::cout << "populate nodes" << std::endl;

  // Loop over cells and add nodes.
  idx_t iNodeLocalOwned = 0;
  idx_t iNodeLocalGhost = nNodesUnique;
  gidx_t iNodeGlobalOwned = 1;
  gidx_t iNodeGlobalGhost = nNodesUnique + 1;


  for (idx_t t = 0; t < nTiles; ++t) {

    // Calc coordinate xy Jacobian wrt ij (this belongs in tile/grid class).
    const auto i0_j0 = static_cast<size_t>(cellLocalIdx(t, 0, 0));
    const auto i1_j0 = static_cast<size_t>(cellLocalIdx(t, 0, 1));
    const auto i0_j1 = static_cast<size_t>(cellLocalIdx(t, 1, 0));

    const auto dx_di = cellXY[i1_j0].x() - cellXY[i0_j0].x();
    const auto dx_dj = cellXY[i0_j1].x() - cellXY[i0_j0].x();
    const auto dy_di = cellXY[i1_j0].y() - cellXY[i0_j0].y();
    const auto dy_dj = cellXY[i0_j1].y() - cellXY[i0_j0].y();

    std:: cout << t << " dx_di " << dx_di << std::endl;
    std:: cout << t << " dx_dj " << dx_dj << std::endl;
    std:: cout << t << " dy_di " << dy_di << std::endl;
    std:: cout << t << " dy_dj " << dy_dj << std::endl << std::endl;

    for (idx_t j = 0; j < N; ++j) {
      for (idx_t i = 0; i < N; ++i) {

        auto addNode = [&](const idx_t iCellLocal,
          const double di, const double dj) {

          // Get cell centre.
          const auto x0 = cellXY[static_cast<size_t>(iCellLocal)].x();
          const auto y0 = cellXY[static_cast<size_t>(iCellLocal)].y();

          // Set node position.
          const auto apparentXY = PointXY{
            x0 + di * dx_di + dj * dx_dj,
            y0 + di * dy_di + dj * dy_dj
          };

          // Use tile class to convert apparent XY to true XY.
          const auto trueXY = gridTiles.tileCubePeriodicity(apparentXY, t);

          //std::cout << apparentXY << std::endl;
          //std::cout << trueXY << std::endl;

          // Node is owned if trueXY = apparentXY. Otherwise, ghost.
          idx_t iNodeLocal = -1;

          if (gridTiles.tileFromXY(trueXY.data()) == t) {

            iNodeLocal = iNodeLocalOwned++;

            mesh::Nodes::Topology::reset(nodeFlags(iNodeLocal));
            nodeGhost(iNodeLocal) = 0;
            nodeGlobalIdx(iNodeLocal) = iNodeGlobalOwned++;

          } else {

            iNodeLocal = iNodeLocalGhost++;

            mesh::Nodes::Topology::set(nodeFlags(iNodeLocal),
              mesh::Nodes::Topology::GHOST);
            nodeGhost(iNodeLocal) = 1;
            nodeGlobalIdx(iNodeLocal) = iNodeGlobalGhost++;
          }

          // Set xy
          nodeXY(iNodeLocal, XX) = apparentXY.x();
          nodeXY(iNodeLocal, YY) = apparentXY.y();


          // Set node lon and lat.
          const PointLonLat lonLat = csGrid.projection().lonlat(trueXY);
          nodeLonLat(iNodeLocal, LON) = lonLat.lon();
          nodeLonLat(iNodeLocal, LAT) = lonLat.lat();


          const auto checkXY = csGrid.projection().xy(lonLat);

          ATLAS_ASSERT(checkXY == trueXY, "lonlat and xy do not match" );

          // Set remote id to local index (all node owned for now).
          nodeRemoteIdx(iNodeLocal) = iNodeLocal;

          // Set partition to zero (for now).
          nodePart(iNodeLocal) = 0;

          // Update local node grid.
          nodeLocalIdx(t, j + static_cast<idx_t>(dj + 0.5),
            i + static_cast<idx_t>(di + 0.5)) = iNodeLocal;
          //std::cout << t << " " << j + static_cast<idx_t>(dj + 0.5) << " " << i + static_cast<idx_t>(di + 0.5) << std::endl;
          //std::cout << std::endl;


          return;
        };


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
        const auto iCellLocal = cellLocalIdx(t, j, i);

        // Every cell sets top right node.
        addNode(iCellLocal, 0.5, 0.5);

        // i == 0 cells set top left node.
        if (i == 0) addNode(iCellLocal, -0.5, 0.5);

        // j == 0 cells set bottom right node.
        if (j == 0) addNode(iCellLocal, 0.5, -0.5);

        // i == 0 and j == 0 cell sets bottom left node.
        if (i == 0 && j == 0) addNode(iCellLocal, -0.5, -0.5);


        // Get list of nodes for cell
        const auto iNodes = std::array<idx_t, 4> {
          nodeLocalIdx(t, j, i),
          nodeLocalIdx(t, j, i + 1),
          nodeLocalIdx(t, j + 1, i + 1),
          nodeLocalIdx(t, j + 1, i)
        };

        // Set node connectivity.
        nodeConnectivity.set(iCellLocal, iNodes.data());

      }
    }
  }

  std::cout << "meshgen done" << std::endl;
  std::cout << iNodeLocalOwned << " " << iNodeLocalGhost << std::endl;
  std::cout << nNodesUnique << " " << nNodesAll << std::endl;


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
