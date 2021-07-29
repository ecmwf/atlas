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
      "with a cubedsphere grid.", Here());
  }

  // Check for proper stagger
  const auto gridName = grid.name();
  const auto gridStagger = gridName.substr(gridName.rfind("-") - 1, 1);


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

  // Get dimensions of grid
  const auto N      = csGrid.N();
  const auto nTiles = csGrid.GetNTiles();

  ATLAS_TRACE("CubedSphereMeshGenerator::generate");
  Log::debug() << "Number of cells per tile edge = " << std::to_string(N) << std::endl;

  // Set tiles.
  // TODO: this needs to be replaced with regular expression matching.
  auto gridTiles = CubedSphereTiles(util::Config("type", "cubedsphere_lfric"));

  // Number of nodes and cells.
  const auto nNodesUnique  = nTiles * N * N + 2;
  const auto nNodesAll     = nTiles * (N + 1) * (N + 1);
  const auto nCells        = nTiles * N * N;


  // Construct mesh cells.
  auto& cells = mesh.cells();
  cells.add(new mesh::temporary::Quadrilateral(), nCells);
  auto cellRemoteIdxArr = array::make_indexview<idx_t, 1>(cells.remote_index());
  auto cellGlobalIdxArr = array::make_view<gidx_t, 1>(cells.global_index());
  auto cellPartArr      = array::make_view<int, 1>(cells.partition());
  auto cellXyArr        = std::vector<PointXY>(static_cast<size_t>(nCells));

  // Set grid shaped arrays of cell local indices.
  array::ArrayT<int> cellLocalIdxData(nTiles, N, N);
  auto cellLocalIdxGrid = array::make_view<idx_t, 3>(cellLocalIdxData);

  // Loop over all cells and set local index
  idx_t cellIdx = 0;
  for (idx_t t = 0; t < nTiles; ++t) {
    for (idx_t j = 0; j < N; ++j) {
      for (idx_t i = 0; i < N; ++i) {
        cellLocalIdxGrid(t, j, i) = cellIdx++;
      }
    }
  }

  // Loop over grid points and set cell global index and xy
  gidx_t gridIdx = 1;
  auto tji = csGrid.tij().begin();
  for (const auto& xy : csGrid.xy()) {

      // Get local index.
      const auto cellLocalIdx =
        cellLocalIdxGrid((*tji).t(), (*tji).j(), (*tji).i());
      ++tji;

      // Set cell-centroid xy.
      cellXyArr[static_cast<size_t>(cellLocalIdx)] = xy;

      // Set remote id to iLocal (all cells owned for now)
      cellRemoteIdxArr(cellLocalIdx) = cellLocalIdx;

      // Set partition using grid global index.
      cellPartArr(cellLocalIdx) = 0 ;

      // Set cell global index to grid global index.
      cellGlobalIdxArr(cellLocalIdx) = gridIdx++;

  }

  // Construct mesh nodes.
  auto& nodes = mesh.nodes();
  nodes.resize(nNodesAll);
  auto nodeRemoteIdxArr = array::make_view<idx_t, 1>(nodes.remote_index());
  auto nodeGlobalIdxArr = array::make_view<gidx_t, 1>(nodes.global_index());
  auto nodePartArr      = array::make_view<int, 1>(nodes.partition());
  auto nodeLocalXyArr   = array::make_view<double, 2>(nodes.xy());
  auto nodeLonLatArr    = array::make_view<double, 2>(nodes.lonlat());
  auto nodeGhostArr     = array::make_view<int, 1>(nodes.ghost());
  auto nodeFlagsArr     = array::make_view<int, 1>(nodes.flags());
  auto nodeRemoteXyArr  = std::vector<PointXY>{};

  // Set grid shaped array of node local indices.
  array::ArrayT<int> nodeLocalIdxData(nTiles, N + 1, N + 1);
  auto nodeLocalIdxGrid = array::make_view<idx_t, 3>(nodeLocalIdxData);

  // Calculate Jacobian of xy wrt ij for each tile so that we can easily
  // switch between the two.
  struct Jacobian {
    double dxByDi{};
    double dxByDj{};
    double dyByDi{};
    double dyByDj{};
    double diByDx{};
    double diByDy{};
    double djByDx{};
    double djByDy{};
    PointXY xy0{};
  };

  idx_t t = 0;
  auto xyJacobians = std::vector<Jacobian>{};
  std::generate_n(std::back_inserter(xyJacobians), nTiles, [&](){

      // Initialise element.
      auto jacElem = Jacobian{};

      // Get cell indices.
      const auto ij00 = static_cast<size_t>(cellLocalIdxGrid(t, 0, 0));
      const auto ij10 = static_cast<size_t>(cellLocalIdxGrid(t, 0, 1));
      const auto ij01 = static_cast<size_t>(cellLocalIdxGrid(t, 1, 0));

      // Calculate Jacobian.
      jacElem.dxByDi = cellXyArr[ij10].x() - cellXyArr[ij00].x();
      jacElem.dxByDj = cellXyArr[ij01].x() - cellXyArr[ij00].x();
      jacElem.dyByDi = cellXyArr[ij10].y() - cellXyArr[ij00].y();
      jacElem.dyByDj = cellXyArr[ij01].y() - cellXyArr[ij00].y();

      // Calculate inverse Jacobian.
      const auto invDet =
        1./(jacElem.dxByDi * jacElem.dyByDj - jacElem.dxByDj * jacElem.dyByDi);
      jacElem.diByDx =  jacElem.dyByDj * invDet;
      jacElem.diByDy = -jacElem.dxByDj * invDet;
      jacElem.djByDx = -jacElem.dyByDi * invDet;
      jacElem.djByDy =  jacElem.dxByDi * invDet;

      // Extrapolate to get xy of node(t, 0, 0)
      jacElem.xy0 = PointXY{
        cellXyArr[ij00].x() - 0.5 * jacElem.dxByDi - 0.5 * jacElem.dxByDj,
        cellXyArr[ij00].y() - 0.5 * jacElem.dyByDi - 0.5 * jacElem.dyByDj
      };

      ++t;
      return jacElem;

    });

  // Initialise node connectivity.
  auto& nodeConnectivity = mesh.cells().node_connectivity();

  // Loop over cells and add nodes.
  idx_t nodeLocalOwnedIdx = 0;
  idx_t nodeLocalGhostIdx = nNodesUnique;
  gidx_t nodeGlobalOwnedIdx = 1;
  gidx_t nodeGlobalGhostIdx = nNodesUnique + 1;

  for (idx_t t = 0; t < nTiles; ++t) {
    for (idx_t jNode = 0; jNode < N + 1; ++jNode) {
      for (idx_t iNode = 0; iNode < N + 1 ; ++iNode) {

        // Get cell indices.
        const auto iCell = std::max(0, iNode - 1);
        const auto jCell = std::max(0, jNode - 1);
        const auto cellLocalIdx = cellLocalIdxGrid(t, jCell, iCell);

        // Get cell centre.
        const auto x0 = cellXyArr[static_cast<size_t>(cellLocalIdx)].x();
        const auto y0 = cellXyArr[static_cast<size_t>(cellLocalIdx)].y();

        // Extrapolate node xy from cell centre.
        const auto di = iNode - iCell - 0.5;
        const auto dj = jNode - jCell - 0.5;
        const auto nodeLocalXy = PointXY{
          x0 + di * xyJacobians[static_cast<size_t>(t)].dxByDi
             + dj * xyJacobians[static_cast<size_t>(t)].dxByDj,
          y0 + di * xyJacobians[static_cast<size_t>(t)].dyByDi
             + dj * xyJacobians[static_cast<size_t>(t)].dyByDj
        };

        // Use tile class to convert local xy to remote xy.
        const auto nodeRemoteXy = gridTiles.tileCubePeriodicity(nodeLocalXy, t);

        // Node is owned if nodeRemoteXy and nodeLocalXy are on same tile.
        idx_t nodeLocalIdx{};
        if (gridTiles.tileFromXY(nodeRemoteXy.data()) == t) {

          // Owned node.

          // Get local node index.
          nodeLocalIdx = nodeLocalOwnedIdx++;

          // Set flags
          mesh::Nodes::Topology::reset(nodeFlagsArr(nodeLocalIdx));
          nodeGhostArr(nodeLocalIdx) = 0;

          // Set global index
          nodeGlobalIdxArr(nodeLocalIdx) = nodeGlobalOwnedIdx++;

          // Set remote id to local index (all node owned for now).
          nodeRemoteIdxArr(nodeLocalIdx) = nodeLocalIdx;

          // Set partition to cell part
          nodePartArr(nodeLocalIdx) = cellPartArr(cellLocalIdx);

        } else {

          // Ghost node.

          // Get local node index.
          nodeLocalIdx = nodeLocalGhostIdx++;

          // Set flags.
          mesh::Nodes::Topology::set(nodeFlagsArr(nodeLocalIdx),
            mesh::Nodes::Topology::GHOST);
          nodeGhostArr(nodeLocalIdx) = 1;

          // Set global index (ghost points have unique global index).
          nodeGlobalIdxArr(nodeLocalIdx) = nodeGlobalGhostIdx++;

          // Need to work this out once we've populated rest of mesh.
          nodeRemoteIdxArr(nodeLocalIdx) = -1;

          // Need to work this out once we've populated rest of mesh.
          nodePartArr(nodeLocalIdx) = -1;

          // Keep track of remote xy.
          nodeRemoteXyArr.push_back(nodeRemoteXy);

        }

        // Set xy
        nodeLocalXyArr(nodeLocalIdx, XX) = nodeLocalXy.x();
        nodeLocalXyArr(nodeLocalIdx, YY) = nodeLocalXy.y();

        // Set node lon and lat.
        const PointLonLat lonLat = csGrid.projection().lonlat(nodeRemoteXy);
        nodeLonLatArr(nodeLocalIdx, LON) = lonLat.lon();
        nodeLonLatArr(nodeLocalIdx, LAT) = lonLat.lat();

        // Update local node grid.
        nodeLocalIdxGrid(t, jNode, iNode) = nodeLocalIdx;

        // Node layout relative to cell.
        //
        //  ^  N3 - N2
        //  ^  |  C  |
        //  j  N0 - N1
        //      i > >
        //
        // C: Cell
        // N0, N1, N2, N3: Nodes
        if (iNode > 0 && jNode > 0) {

          // Get nodes of quadrilateral cell.
          const auto quadNodes = std::array<idx_t, 4> {
            nodeLocalIdxGrid(t, jNode - 1, iNode - 1),
            nodeLocalIdxGrid(t, jNode - 1, iNode    ),
            nodeLocalIdxGrid(t, jNode    , iNode    ),
            nodeLocalIdxGrid(t, jNode    , iNode - 1)
          };

          // Set node connectivity.
          nodeConnectivity.set(cellLocalIdx, quadNodes.data());
        }

      }
    }
  }

  // Set remote indices of ghost points.
  auto nodeRemoteXyIt = nodeRemoteXyArr.begin();
  for (idx_t nodeLocalIdx = nNodesUnique;
    nodeLocalIdx < nNodesAll; ++nodeLocalIdx) {

    // Get remote xy
    const auto nodeRemoteXy = *nodeRemoteXyIt++;

    // Get remote t
    const auto t = static_cast<size_t>(
      gridTiles.tileFromXY(nodeRemoteXy.data()));

    // Get remote i and j
    const auto dx = nodeRemoteXy.x() - xyJacobians[t].xy0.x();
    const auto dy = nodeRemoteXy.y() - xyJacobians[t].xy0.y();

    // Round to deal with potential floating point error.
    const auto i = static_cast<idx_t>(
      std::round(dx * xyJacobians[t].diByDx + dy * xyJacobians[t].diByDy));
    const auto j = static_cast<idx_t>(
      std::round(dx * xyJacobians[t].djByDx + dy * xyJacobians[t].djByDy));

    // Set remote index and partition.
    const auto nodeRemoteIdx = nodeLocalIdxGrid(t, j, i);
    nodeRemoteIdxArr(nodeLocalIdx) = nodeRemoteIdx;
    nodePartArr(nodeLocalIdx) = nodePartArr(nodeRemoteIdx);

  }

  nodes.metadata().set( "parallel", true );

  // check that node counts are correct
  ATLAS_ASSERT(nodeLocalOwnedIdx = nNodesUnique);
  ATLAS_ASSERT(nodeLocalGhostIdx = nNodesAll);

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
