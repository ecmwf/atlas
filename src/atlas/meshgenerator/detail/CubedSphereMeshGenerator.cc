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

  // Check for correct grid and need for mesh
  ATLAS_ASSERT(!mesh.generated());
  if (!CubedSphereGrid(grid)) {
    throw_Exception("CubedSphereMeshGenerator can only work "
    "with a cubedsphere grid.", Here());
  }

  // Check for correct stagger
  const auto gridName = grid.name();
  const auto gridStagger = gridName.substr(gridName.rfind("-") - 1, 1);


  if (gridStagger != "C") {
    throw_Exception("CubedSphereMeshGenerator can only work with a"
    "cell-centroid grid. Try FV3CubedSphereMeshGenerator instead.");
  }

  // Partitioner
  const auto partitioner = grid::Partitioner("equal_regions", util::Config("coordinates", "lonlat"));
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

  const auto nNodesUnique = nTiles * N * N + 2;
  const auto nNodesAll    = nTiles * (N + 1) * (N + 1);
  const auto nCells       = nTiles * N * N;

  // Define bad index values.
  constexpr auto badIdx = std::numeric_limits<idx_t>::min();
  constexpr auto badGlobalIdx = std::numeric_limits<gidx_t>::min();

  ATLAS_TRACE("CubedSphereMeshGenerator::generate");
  Log::debug() << "Number of cells per tile edge = "
    << std::to_string(N) << std::endl;

  // Set tiles.
  // TODO: this needs to be replaced with regular expression matching.
  auto gridTiles = CubedSphereTiles(util::Config("type", "cubedsphere_lfric"));

  // Get partition information.
  const auto nParts =   mpi::comm().size();
  const auto thisPart = mpi::comm().rank();

  // Helper functions to get node and cell idx from (t, j, i).
  const auto getNodeIdx = [&](idx_t t, idx_t j, idx_t i){
    // Adjust bounds.
    t = std::max(std::min(t, nTiles - 1), 0);
    j = std::max(std::min(j, N), 0);
    i = std::max(std::min(i, N), 0);
    return static_cast<size_t>(t * (N + 1) * (N + 1) + j * (N + 1) + i);
  };

  const auto getCellIdx = [&](idx_t t, idx_t j, idx_t i){
    // Adjust bounds.
    t = std::max(std::min(t, nTiles - 1), 0);
    j = std::max(std::min(j, N - 1), 0);
    i = std::max(std::min(i, N - 1), 0);
    return static_cast<size_t>(t * N * N + j * N + i);
  };

  // Define cell record.
  struct CellRecord {
    gidx_t  globalIdx{badGlobalIdx};
    idx_t   remoteIdx{badIdx};
    idx_t   part{badIdx};
    PointXY xy{};
  };

  // ij bounding box for each face (this partition only).
  struct BoundingBox {
    idx_t iBegin{std::numeric_limits<idx_t>::max()};
    idx_t iEnd{std::numeric_limits<idx_t>::min()};
    idx_t jBegin{std::numeric_limits<idx_t>::max()};
    idx_t jEnd{std::numeric_limits<idx_t>::min()};
  };

  // Make list of all cells.
  auto globalCells = std::vector<CellRecord>(static_cast<size_t>(nCells));

  // Initialise bounding box.
  auto cellBounds = std::vector<BoundingBox>(static_cast<size_t>(nTiles));

  // Loop over grid.
  auto tjiIt = csGrid.tij().begin();
  auto xyIt = csGrid.xy().begin();
  auto cellRemoteIdx = std::vector<idx_t>(static_cast<size_t>(nParts));
  idx_t nCellsLocal = 0;

  std::cout << "Coping grid to cells" << std::endl;

  for (gidx_t gridIdx = 1; gridIdx < csGrid.size() + 1; ++gridIdx) {

    // Get cell index
    const auto t = (*tjiIt).t();
    const auto j = (*tjiIt).j();
    const auto i = (*tjiIt).i();
    const auto cellIdx = getCellIdx(t, j, i);
    auto& cell = globalCells[cellIdx];

    // Set global index.
    cell.globalIdx = gridIdx;

    // Set partition and remote index.
    cell.part = distribution.partition(gridIdx - 1);

    //cell.part = static_cast<idx_t>(thisPart);

    cell.remoteIdx = cellRemoteIdx[static_cast<size_t>(cell.part)]++;

    // Set xy
    cell.xy = *xyIt;

    if (cell.part == static_cast<idx_t>(thisPart)) {

      // Keep track of local (t, j, i) bounds.
      auto& bounds = cellBounds[static_cast<size_t>(t)];
      bounds.iBegin = std::min(bounds.iBegin, i    );
      bounds.iEnd   = std::max(bounds.iEnd  , i + 1);
      bounds.jBegin = std::min(bounds.jBegin, j    );
      bounds.jEnd   = std::max(bounds.jEnd  , j + 1);

      // Count number of local cells.
      ++nCellsLocal;

    }
    // Increment iterators.
    ++tjiIt;
    ++xyIt;
  }

  std::cout << "Setting jacobian" << std::endl;

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
    PointXY xy00{};
  };

  idx_t t = 0;
  auto xyJacobians = std::vector<Jacobian>{};
  std::generate_n(std::back_inserter(xyJacobians), nTiles, [&](){

      // Initialise element.
      auto jacElem = Jacobian{};

      // Get cell positions.
      const auto& xy00 = globalCells[getCellIdx(t, 0, 0)].xy;
      const auto& xy10 = globalCells[getCellIdx(t, 0, 1)].xy;
      const auto& xy01 = globalCells[getCellIdx(t, 1, 0)].xy;

      // Calculate Jacobian.
      jacElem.dxByDi = xy10.x() - xy00.x();
      jacElem.dxByDj = xy01.x() - xy00.x();
      jacElem.dyByDi = xy10.y() - xy00.y();
      jacElem.dyByDj = xy01.y() - xy00.y();

      // Calculate inverse Jacobian.
      const auto invDet =
        1./(jacElem.dxByDi * jacElem.dyByDj - jacElem.dxByDj * jacElem.dyByDi);
      jacElem.diByDx =  jacElem.dyByDj * invDet;
      jacElem.diByDy = -jacElem.dxByDj * invDet;
      jacElem.djByDx = -jacElem.dyByDi * invDet;
      jacElem.djByDy =  jacElem.dxByDi * invDet;

      // Extrapolate cell(t, 0, 0) xy to get node(t, 0, 0) xy.
      jacElem.xy00 = PointXY{
        xy00.x() - 0.5 * jacElem.dxByDi - 0.5 * jacElem.dxByDj,
        xy00.y() - 0.5 * jacElem.dyByDi - 0.5 * jacElem.dyByDj
      };

      ++t;
      return jacElem;

    });

  // Define node record.
  struct NodeRecord {
    gidx_t  globalIdx{badGlobalIdx};
    idx_t   remoteIdx{badIdx};
    idx_t   localIdx{badIdx};
    idx_t   remotePart{badIdx};
    idx_t   localPart{badIdx};
    idx_t   t{badIdx};
    PointXY localXy{};
    PointXY remoteXy{};
  };


  std::cout << "Building global nodes" << std::endl;

  auto globalNodes = std::vector<NodeRecord>(static_cast<size_t>(nNodesAll));

  // Loop over *all* nodes.
  auto nNodesOwned = std::vector<idx_t>(nParts);
  gidx_t nodeGlobalOwnedIdx = 1;
  gidx_t nodeGlobalGhostIdx = nNodesUnique + 1;

  auto edgeGhostNodes = std::vector<NodeRecord*>{};

  for (idx_t t = 0; t < nTiles; ++t) {
    for (idx_t jNode = 0; jNode < N + 1; ++jNode) {
      for (idx_t iNode = 0; iNode < N + 1; ++iNode) {

        // Set and get node.
        auto& node = globalNodes[getNodeIdx(t, jNode, iNode)];

        // Get owning cell.
        const auto iCell = std::max(iNode - 1, 0);
        const auto jCell = std::max(jNode - 1, 0);
        const auto& cell = globalCells[getCellIdx(t, jCell, iCell)];

        // Extrapolate local xy from cell centre.
        const auto di = iNode - iCell - 0.5;
        const auto dj = jNode - jCell - 0.5;
        const auto& jac = xyJacobians[static_cast<size_t>(t)];
        node.localXy = PointXY{
            cell.xy.x() + di * jac.dxByDi + dj * jac.dxByDj,
            cell.xy.y() + di * jac.dyByDi + dj * jac.dyByDj
          };

        // Get remote xy from tile class and "implied" t.
        node.remoteXy = gridTiles.tileCubePeriodicity(node.localXy, t);

        // Local part always taken from owning cell.
        node.localPart = cell.part;

        // Get actual t from tile class.
        node.t = gridTiles.tileFromXY(node.remoteXy.data());

        // Node is an edge ghost point if actual and implied t differ.
        if (t == node.t) {

          // Non edge ghost.

          // Set global index.
          node.globalIdx = nodeGlobalOwnedIdx++;

          // Set parition and remote index.
          node.remotePart = node.localPart;
          node.remoteIdx = nNodesOwned[static_cast<size_t>(node.remotePart)]++;

        } else {

          // Edge ghost.


          // Set global index
          node.globalIdx = nodeGlobalGhostIdx++;


          //Add to list and sort later.
          edgeGhostNodes.push_back(&node);

        }
      }
    }
  }

  std::cout << "sort edge-ghosts" << std::endl;

  // Sort out edge-ghost nodes.
  for (auto nodePtr : edgeGhostNodes) {

    // Get jacobian.
    std::cout << "t " << nodePtr->t << std::endl;
    const auto& jac = xyJacobians[static_cast<size_t>(nodePtr->t)];

    // Get actual i and j.
    const auto dx = nodePtr->remoteXy.x() - jac.xy00.x();
    const auto dy = nodePtr->remoteXy.y() - jac.xy00.y();
    const auto i = static_cast<idx_t>(
      std::round(dx * jac.diByDx + dy * jac.diByDy));
    const auto j = static_cast<idx_t>(
      std::round(dx * jac.djByDx + dy * jac.djByDy));

    // Get remote node.

    std::cout << t << " " << j << " " << i << std::endl;
    const auto& remoteNode = globalNodes[getNodeIdx(nodePtr->t, j, i)];

    // Set partition and remote index.
    nodePtr->remotePart = remoteNode.localPart;
    nodePtr->remoteIdx = remoteNode.globalIdx - 1;

  }



  // Make list of local owned and ghost nodes.
  auto localNodesOwned = std::vector<const NodeRecord*>{};
  auto localNodesGhost = std::vector<const NodeRecord*>{};

  // Loop over *local* nodes.
  idx_t localNodeOwnedIdx = 0;
  idx_t localNodeGhostIdx = nNodesOwned[thisPart];


  std::cout << "Selecting local nodes" << std::endl;

  for (idx_t t = 0; t < nTiles; ++t) {

    const auto bounds = cellBounds[static_cast<size_t>(t)];

    for (idx_t jNode = 0; jNode < N + 1; ++jNode) {
      for (idx_t iNode = 0; iNode < N + 1; ++iNode) {

        // Get node.
        const auto nodeIdx = getNodeIdx(t, jNode, iNode);
        auto& node = globalNodes[nodeIdx];

        // Check if node is owned or edge ghost point.
        if (node.localPart == static_cast<idx_t>(thisPart)) {

          if (node.t == t) {

            // Node is owned.

            // Set local index.
            node.localIdx = localNodeOwnedIdx++;

            std::cout << node.remoteIdx << " " << node.localIdx << std::endl;

            // Add node to list.
            localNodesOwned.push_back(&node);

          } else {

            // Node is an edge ghost point.

            // Set local index.
            node.localIdx = localNodeGhostIdx++;

            node.remoteIdx = node.globalIdx - 1;


            // Add node to list.
            localNodesGhost.push_back(&node);

          }
        } else {

          // Might be a ghost point. Need to check non-owning surrounding cells.
          const auto& cell1 = globalCells[getCellIdx(t, jNode    , iNode - 1)];
          const auto& cell2 = globalCells[getCellIdx(t, jNode    , iNode    )];
          const auto& cell3 = globalCells[getCellIdx(t, jNode - 1, iNode    )];

          if ( cell1.part == static_cast<idx_t>(thisPart) ||
               cell2.part == static_cast<idx_t>(thisPart) ||
               cell3.part == static_cast<idx_t>(thisPart)) {

            // Set local index.
            node.localIdx = localNodeGhostIdx++;

            // Add node to list.
            localNodesGhost.push_back(&node);

          }
        }
      }
    }
  }

  std::cout << "writing to nodes to mesh." << std::endl;

  // We now have enough information to construct mesh.
  const auto nNodesLocalOwned = static_cast<idx_t>(localNodesOwned.size());
  const auto nNodesLocalGhost = static_cast<idx_t>(localNodesGhost.size());
  const auto nNodesLocalAll = nNodesLocalOwned + nNodesLocalGhost;

  // Resize nodes.
  mesh.nodes().resize(nNodesLocalAll);

  // Get field views
  auto meshNodesXy        = array::make_view<double, 2>(mesh.nodes().xy());
  auto meshNodesLonLat    = array::make_view<double, 2>(mesh.nodes().lonlat());
  auto meshNodesGobalIdx  = array::make_view<gidx_t, 1>(mesh.nodes().global_index());
  auto meshNodesRemoteIdx = array::make_indexview<idx_t, 1>(mesh.nodes().remote_index());
  auto meshNodesPart      = array::make_view<int, 1>(mesh.nodes().partition());
  auto meshNodesGhost     = array::make_view<int, 1>(mesh.nodes().ghost());
  auto meshNodesFlags     = array::make_view<int, 1>(mesh.nodes().flags());

  // Set owned nodes.
  auto nodeOwnedIt = localNodesOwned.begin();
  auto nodeGhostIt = localNodesGhost.begin();
  for (idx_t nodeIdx = 0; nodeIdx < nNodesLocalAll; ++nodeIdx) {

    // Get node record.
    const auto ghost = nodeIdx >= nNodesLocalOwned;
    const auto& node = ghost ? **nodeGhostIt++ : **nodeOwnedIt++;

    // Set xy.
    meshNodesXy(nodeIdx, XX) = node.localXy.x();
    meshNodesXy(nodeIdx, YY) = node.localXy.y();

    // Set lon-lat.
    const auto lonLat = csGrid.projection().lonlat(node.remoteXy);
    meshNodesLonLat(nodeIdx, LON) = lonLat.lon();
    meshNodesLonLat(nodeIdx, LAT) = lonLat.lat();

    // Set global index.
    meshNodesGobalIdx(nodeIdx) = node.globalIdx;

    // Set remote index
    meshNodesRemoteIdx(nodeIdx) = node.remoteIdx;
    if (!ghost) meshNodesRemoteIdx(nodeIdx) = - 1;


    // Set partition.
    meshNodesPart(nodeIdx) = node.remotePart;

    // Set ghost flag.
    meshNodesGhost(nodeIdx) = ghost;

    mesh::Nodes::Topology::reset(meshNodesFlags(nodeIdx));
    if (ghost) mesh::Nodes::Topology::set(meshNodesFlags(nodeIdx), mesh::Nodes::Topology::GHOST);
  }

  std::cout << "writing to cells to mesh. " << nCellsLocal << std::endl;

  // Resize cells.
  mesh.cells().add(new mesh::temporary::Quadrilateral(), nCellsLocal);

  // Set field views.
  auto meshCellsRemoteIdx = array::make_indexview<idx_t, 1>(mesh.cells().remote_index());
  auto meshCellsGlobalIdx = array::make_view<gidx_t, 1>(mesh.cells().global_index());
  auto meshCellsPart      = array::make_view<int, 1>(mesh.cells().partition());

  // Set local cells.
  auto& nodeConnectivity = mesh.cells().node_connectivity();
  idx_t cellIdx = 0;

  for (idx_t t = 0; t < nTiles; ++t) {

    const auto bounds = cellBounds[static_cast<size_t>(t)];

    for (idx_t jCell = 0; jCell < N; ++jCell) {
      for (idx_t iCell = 0; iCell < N; ++iCell) {

        // Get cell.
        const auto& cell = globalCells[getCellIdx(t, jCell, iCell)];

        // Only add cells on this partition.
        if (cell.part == static_cast<idx_t>(thisPart)) {

          // Set global index.
          meshCellsGlobalIdx(cellIdx) = cell.globalIdx;

          // Set remote index.
          meshCellsRemoteIdx(cellIdx) = cell.remoteIdx;

          // Set partition.
          meshCellsPart(cellIdx) = cell.part;

          // Set quadrilateral.
          const auto i0 = globalNodes[getNodeIdx(t, jCell    , iCell    )].localIdx;
          const auto i1 = globalNodes[getNodeIdx(t, jCell    , iCell + 1)].localIdx;
          const auto i2 = globalNodes[getNodeIdx(t, jCell + 1, iCell + 1)].localIdx;
          const auto i3 = globalNodes[getNodeIdx(t, jCell + 1, iCell    )].localIdx;
          const auto quadNodeIdx = std::array<idx_t, 4> {i0, i1, i2, i3};


          std::cout << quadNodeIdx << std::endl;

          // Set connectivity.
          nodeConnectivity.set(cellIdx, quadNodeIdx.data());

          // Incriment cell index.
          ++cellIdx;

        }
      }
    }
  }

  //generateGlobalElementNumbering( mesh );
  mesh.nodes().metadata().set( "periodic", true );
  mesh.nodes().metadata().set( "parallel", true );

  std::cout << "mesh generated " << cellIdx << " " <<nCellsLocal <<std::endl;

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
