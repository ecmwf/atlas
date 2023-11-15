/*
 * (C) Crown Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <string>
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
#include "atlas/meshgenerator/detail/cubedsphere/CubedSphereUtility.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/projection/detail/CubedSphereProjectionBase.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"
#include "atlas/util/Constants.h"
#include "atlas/util/CoordinateEnums.h"

#define DEBUG_OUTPUT 0
#define DEBUG_OUTPUT_DETAIL 0

namespace atlas {
namespace meshgenerator {

// -----------------------------------------------------------------------------

CubedSphereMeshGenerator::CubedSphereMeshGenerator(const eckit::Parametrisation& p) {
    // Get mpi_comm
    std::string mpi_comm = mpi::comm().name();
    p.get("mpi_comm", mpi_comm);
    options.set("mpi_comm", mpi_comm);

    configure_defaults();

    // Get number of partitions.
    int nb_parts;
    if (p.get("nb_parts", nb_parts)) {
        options.set("nb_parts", nb_parts);
    }

    // Get this partition.
    int part;
    if (p.get("part", part)) {
        options.set("part", part);
    }

    // Get halo size.
    int halo;
    if (p.get("halo", halo)) {
        options.set("halo", halo);
    }

    // Get partitioner.
    std::string partitioner;
    if (p.get("partitioner", partitioner) && partitioner.size()) {
        options.set("partitioner", partitioner);
    }
}

// -----------------------------------------------------------------------------


void CubedSphereMeshGenerator::configure_defaults() {
    auto& comm = mpi::comm(options.getString("mpi_comm"));

    // This option sets number of partitions.
    options.set("nb_parts", comm.size());

    // This option sets the part that will be generated.
    options.set("part", comm.rank());

    // This options sets the number of halo elements around each partition.
    options.set("halo", 0);

    // This options sets the default partitioner.
    options.set<std::string>("partitioner", "cubedsphere");
}

// -----------------------------------------------------------------------------

void CubedSphereMeshGenerator::generate(const Grid& grid, Mesh& mesh) const {
    // Get partitioner type and number of partitions from config.
    const idx_t nParts         = static_cast<idx_t>(options.get<size_t>("nb_parts"));
    const std::string partType = options.get<std::string>("partitioner");

    auto partConfig = util::Config{};
    partConfig.set("type", partType);
    partConfig.set("partitions", nParts);

    // Use lonlat instead of xy for non cubedsphere partitioner.
    if (partType != "cubedsphere") {
        partConfig.set("coordinates", "lonlat");
    }

    // Set distribution.
    mpi::Scope mpi_scope(options.getString("mpi_comm"));
    const auto partitioner  = grid::Partitioner(partConfig);
    const auto distribution = grid::Distribution(grid, partitioner);

    generate(grid, distribution, mesh);
}

// -----------------------------------------------------------------------------

void CubedSphereMeshGenerator::generate(const Grid& grid, const grid::Distribution& distribution, Mesh& mesh) const {
    // Check for correct grid and need for mesh
    ATLAS_ASSERT(!mesh.generated());

    // Cast grid to cubed sphere grid.
    const auto csGrid = CubedSphereGrid(grid);

    // Check for successful cast.
    if (!csGrid) {
        throw_Exception(
            "CubedSphereMeshGenerator can only work "
            "with a cubedsphere grid.",
            Here());
    }

    // Check for correct grid stagger.
    if (csGrid.stagger() != "C") {
        throw_Exception(
            "CubedSphereMeshGenerator will only work with a "
            "cell-centroid grid.",
            Here());
    }

    // Check for sensible halo size.
    if (options.get<idx_t>("halo") > csGrid.N()) {
        throw_Exception("Halo size " + std::to_string(options.get<idx_t>("halo")) +
                            " "
                            "is larger than grid size " +
                            std::to_string(csGrid.N()) + ".",
                        Here());
    }

    // Clone some grid properties.
    setGrid(mesh, csGrid, distribution);
    mesh.metadata().set("mpi_comm",options.getString("mpi_comm"));

    generate_mesh(csGrid, distribution, mesh);
}

// -----------------------------------------------------------------------------

void CubedSphereMeshGenerator::generate_mesh(const CubedSphereGrid& csGrid, const grid::Distribution& distribution,
                                             Mesh& mesh) const {
    ATLAS_TRACE("CubedSphereMeshGenerator::generate");

    // ---------------------------------------------------------------------------
    // CUBED SPHERE MESH GENERATOR
    // ---------------------------------------------------------------------------
    //
    // Mesh generator creates a cubed sphere mesh by generating individual meshes
    // for each cubed sphere tile and then working out the correspondence between
    // overlapping nodes and cells.
    //
    // Meshgenerator places cell at each grid point on this partition. Halo
    // cells are added to the mesh if options.get("halo") > 0. The halo cells may
    // either be interior to the tile, || exterior. Interior halo cells share
    // their xy coordinate and global ID with the corresponding cells on other
    // partitions. Exterior halo cells have a unique global ID and their xy
    // coordinates are extrapolated from the interior tile. The global IDs of non-
    // halo cells match the global IDs of the cell-centre grid points.
    //
    // Nodes are added around cells. The nodes around interior cells are assigned
    // an owner by the following rules:
    //    * node (i > 0, j > 0) is owned by cell (i - 1, j - 1)
    //    * node (i = 0, j > 0) is owned by cell (0    , j - 1)
    //    * node (i > 0, j = 0) is owned by cell (i - 1, 0    )
    //    * node (i = 0, j = 0) is owned by cell (0    , 0    )
    //
    // The partition of the owning cell determines the partition of the node.
    // Ghost nodes are added to the mesh to complete cells at partition
    // boundaries, cells exterior to the tiles, or on tile edges which do not
    // own nodes.
    //
    // There are several stages to the mesh generator:
    //    1. Preamble.
    //    2. Define global cell distribution.
    //    3. Define global node distribution.
    //    4. Locate local cells.
    //    5. Locate local nodes
    //    6. Assign nodes to mesh.
    //    7. Assign cells to mesh.
    //    8. Finalise.

    // ---------------------------------------------------------------------------
    // 1. PREAMBLE
    //    Setup some general parameters of the mesh, as well as define useful
    //    lambda functions and structs.
    // ---------------------------------------------------------------------------

    using Topology = atlas::mesh::Nodes::Topology;
    using atlas::array::make_datatype;
    using atlas::array::make_shape;

    using namespace detail::cubedsphere;

    struct GlobalElem;
    struct LocalElem;

    // Define bad index values.
    constexpr idx_t undefinedIdx       = -1;
    constexpr idx_t undefinedGlobalIdx = -1;

    // Define cell/node type.
    enum class ElemType
    {
        UNDEFINED,  // Not set. Writing this type to mesh is an error.
        OWNER,      // Owner element on this partition.
        HALO        // Ghost or halo element.
    };

    // Struct to store the information of globals cells/nodes.
    // These must be generated for all cells/nodes on all partitions. Data is
    // stored in (t, j, i) row-major order.
    struct GlobalElem {
        LocalElem* localPtr{};                 // Pointer to local element.
        gidx_t globalIdx{undefinedGlobalIdx};  // Global index.
        idx_t remoteIdx{undefinedIdx};         // Remote index.
        int part{undefinedIdx};                // Partition.
    };

    // Struct to store additional information of local cells/nodes.
    // These are only generated for nodes on this partition. Data is stored as a
    // list which can be sorted.
    struct LocalElem {
        GlobalElem* globalPtr{};             // Pointer to global element.
        ElemType type{ElemType::UNDEFINED};  // Cell/node Type.
        int halo{undefinedIdx};              // Halo level.
        idx_t t{};                           // t, i and j.
        idx_t i{};
        idx_t j{};
    };

    // Define ij bounding box for each face (this partition only).
    struct BoundingBox {
        idx_t iBegin{std::numeric_limits<idx_t>::max()};
        idx_t iEnd{std::numeric_limits<idx_t>::min()};
        idx_t jBegin{std::numeric_limits<idx_t>::max()};
        idx_t jEnd{std::numeric_limits<idx_t>::min()};
    };

    // Get dimensions of grid
    const idx_t N = csGrid.N();

    // Get size of halo.
    int nHalo = 0;
    options.get("halo", nHalo);

    // Unique non-halo nodes and cells.
    const idx_t nNodesUnique = 6 * N * N + 2;
    const idx_t nCellsUnique = 6 * N * N;

    // Total array sizes (including invalid corner ijs).
    const idx_t nNodesArray = 6 * (N + 2 * nHalo + 1) * (N + 2 * nHalo + 1);
    const idx_t nCellsArray = 6 * (N + 2 * nHalo) * (N + 2 * nHalo);

    // Total number of possible cells and nodes (including halos).
    const idx_t nNodesTotal = nNodesArray - 6 * 4 * nHalo * nHalo;
    const idx_t nCellsTotal = nCellsArray - 6 * 4 * nHalo * nHalo;

    // Projection and jacobian.
    const auto& csProjection = csGrid.cubedSphereProjection();
    const auto jacobian      = NeighbourJacobian(csGrid);

    // Get partition information.
    const int nParts   = options.get<int>("nb_parts");
    const int thisPart = options.get<int>("part");

    // Define an index counter.
    const auto idxSum = [](const std::vector<idx_t>& idxCounts) -> idx_t {
        return std::accumulate(idxCounts.begin(), idxCounts.end(), 0);
    };

    // Helper functions to get node vector idx from (i, j, t).
    const auto getNodeIdx = [&](idx_t i, idx_t j, idx_t t) -> size_t {
        const idx_t rowSize  = (N + 2 * nHalo + 1);
        const idx_t tileSize = rowSize * rowSize;

        return static_cast<size_t>(t * tileSize + (j + nHalo) * rowSize + i + nHalo);
    };

    // Helper functions to get cell vector idx from (i, j, t).
    const auto getCellIdx = [&](idx_t i, idx_t j, idx_t t) -> size_t {
        const idx_t rowSize  = (N + 2 * nHalo);
        const idx_t tileSize = rowSize * rowSize;

        return static_cast<size_t>(t * tileSize + (j + nHalo) * rowSize + i + nHalo);
    };

    // Return true for nodes interior to tile (excluding edge).
    const auto interiorNode = [&](idx_t i, idx_t j) -> bool { return i > 0 && i < N && j > 0 && j < N; };

    // return true for nodes exterior to tile (excluding edge).
    const auto exteriorNode = [&](idx_t i, idx_t j) -> bool { return i < 0 || i > N || j < 0 || j > N; };

    // Return true for cells interior to tile.
    const auto interiorCell = [&](idx_t i, idx_t j) -> bool { return i >= 0 && i < N && j >= 0 && j < N; };

    // Return true for impossible node (i, j) combinations.
    const auto invalidNode = [&](idx_t i, idx_t j) -> bool {
        const bool inCorner = (i < 0 && j < 0) ||  // Bottom-left corner.
                              (i > N && j < 0) ||  // Bottom-right corner.
                              (i > N && j > N) ||  // Top-right corner.
                              (i < 0 && j > N);    // Top-left corner.

        if (inCorner) {
            return true;
        }

        const bool outOfBounds = i < -nHalo || i > N + nHalo || j < -nHalo || j > N + nHalo;

        if (outOfBounds) {
            return true;
        }

        return false;
    };

    // Return true for impossible cell (i, j) combinations.
    const auto invalidCell = [&](idx_t i, idx_t j) -> bool {
        const bool inCorner = (i < 0 && j < 0) ||          // Bottom-left corner.
                              (i > N - 1 && j < 0) ||      // Bottom-right corner.
                              (i > N - 1 && j > N - 1) ||  // Top-right corner.
                              (i < 0 && j > N - 1);        // Top-left corner.

        if (inCorner) {
            return true;
        }

        const bool outOfBounds = i < -nHalo || i > N + nHalo - 1 || j < -nHalo || j > N + nHalo - 1;

        if (outOfBounds) {
            return true;
        }

        return false;
    };

    // Return the (i, j) of cell that owns node (i, j).
    // It's possible that this may need to be user-defined in the future.
    const auto nodeOwnerCell = [](idx_t iNode, idx_t jNode) -> std::pair<idx_t, idx_t> {
        return std::make_pair(iNode > 0 ? iNode - 1 : iNode, jNode > 0 ? jNode - 1 : jNode);
    };

    // ---------------------------------------------------------------------------
    // 2. GLOBAL CELL DISTRIBUTION
    //    Need to construct a lightweight global mesh. This will help us pair up
    //    local and remote indices on different partitions.
    // ---------------------------------------------------------------------------

    // Make vector of all possible cells.
    auto globalCells = std::vector<GlobalElem>(static_cast<size_t>(nCellsArray));

    // Initialise bounding box.
    auto cellBounds = std::vector<BoundingBox>(6);

    // Set tij grid iterator.
    auto tijIt = csGrid.tij().begin();

    for (gidx_t gridIdx = 0; gridIdx < csGrid.size(); ++gridIdx) {
        // It's more than likely that the grid order follows the same (t, j, i) row-
        // major order of the mesh. However, we shouldn't assume this is true.

        // Get cell index.
        const idx_t t          = (*tijIt).t();
        const idx_t i          = (*tijIt).i();
        const idx_t j          = (*tijIt).j();
        const size_t cellIdx   = getCellIdx(i, j, t);
        GlobalElem& globalCell = globalCells[cellIdx];

        // Set global index.
        globalCell.globalIdx = gridIdx + 1;

        // Set partition.
        globalCell.part = distribution.partition(gridIdx);

        if (globalCell.part == static_cast<idx_t>(thisPart)) {
            // Keep track of local (t, j, i) bounds.
            BoundingBox& bounds = cellBounds[static_cast<size_t>(t)];
            bounds.iBegin       = std::min(bounds.iBegin, i - nHalo);
            bounds.iEnd         = std::max(bounds.iEnd, i + nHalo + 1);
            bounds.jBegin       = std::min(bounds.jBegin, j - nHalo);
            bounds.jEnd         = std::max(bounds.jEnd, j + nHalo + 1);
        }
        // Increment tij.
        ++tijIt;
    }

    // Set counters for cell local indices.
    auto cellLocalIdxCount = std::vector<idx_t>(nParts, 0);

    // Give possible edge-halo cells a unique global ID.
    gidx_t cellGlobalIdxCount = nCellsUnique + 1;

    for (idx_t t = 0; t < 6; ++t) {
        for (idx_t j = -nHalo; j < N + nHalo; ++j) {
            for (idx_t i = -nHalo; i < N + nHalo; ++i) {
                // Skip invalid cell.
                if (invalidCell(i, j)) {
                    continue;
                }

                // Get cell.
                GlobalElem& globalCell = globalCells[getCellIdx(i, j, t)];

                // Set cell remote index if interior cell.
                if (interiorCell(i, j)) {
                    // Set remote index.
                    globalCell.remoteIdx = cellLocalIdxCount[static_cast<size_t>(globalCell.part)]++;
                }
                else {
                    // Set global index.
                    globalCell.globalIdx = cellGlobalIdxCount++;
                    // Leave remote index and part undefined.
                }
            }
        }
    }

    ATLAS_ASSERT(idxSum(cellLocalIdxCount) == nCellsUnique);
    ATLAS_ASSERT(cellGlobalIdxCount == nCellsTotal + 1);

    // ---------------------------------------------------------------------------
    // 3. GLOBAL NODE DISTRIBUTION
    //    Construct a lightweight global distribution of nodes. Again, this will
    //    help us match up local and remote indices accross partitions.
    // ---------------------------------------------------------------------------

    // Make list of all nodes.
    auto globalNodes = std::vector<GlobalElem>(static_cast<size_t>(nNodesArray));

    // Set counters for local node indices.
    auto nodeLocalIdxCount = std::vector<idx_t>(nParts, 0);

    // Set counter for global indices.
    gidx_t nodeGlobalOwnedIdxCount = 1;
    idx_t nodeGlobalGhostIdxCount  = nNodesUnique + 1;

    for (idx_t t = 0; t < 6; ++t) {
        for (idx_t j = -nHalo; j < N + nHalo + 1; ++j) {
            for (idx_t i = -nHalo; i < N + nHalo + 1; ++i) {
                // Skip if not a valid node.
                if (invalidNode(i, j)) {
                    continue;
                }

                // Get this node.
                GlobalElem& globalNode = globalNodes[getNodeIdx(i, j, t)];

                // Get owner cell.
                idx_t iCell, jCell;
                std::tie(iCell, jCell)      = nodeOwnerCell(i, j);
                const GlobalElem& ownerCell = globalCells[getCellIdx(iCell, jCell, t)];

                if (interiorNode(i, j)) {
                    // Node is definitely an owner.
                    globalNode.globalIdx = nodeGlobalOwnedIdxCount++;
                    globalNode.part      = ownerCell.part;
                    globalNode.remoteIdx = nodeLocalIdxCount[static_cast<size_t>(globalNode.part)]++;
                }
                else if (exteriorNode(i, j)) {
                    // Node is definitely a ghost.
                    globalNode.globalIdx = nodeGlobalGhostIdxCount++;
                    // Leave remote index and partition undefined.
                }
                else {
                    // We're not sure (i.e., node is on a tile edge).

                    // Check that xy is on this tile.
                    PointXY xy = jacobian.xy(PointIJ(i, j), t);
                    xy         = jacobian.snapToEdge(xy, t);

                    // This will only determine if tGlobal does not match t.
                    // This is cheaper than determining the correct tGlobal.
                    idx_t tGlobal = csProjection.getCubedSphereTiles().indexFromXY(xy.data());

                    if (tGlobal == t) {
                        // Node is an owner.
                        globalNode.globalIdx = nodeGlobalOwnedIdxCount++;
                        globalNode.part      = ownerCell.part;
                        globalNode.remoteIdx = nodeLocalIdxCount[static_cast<size_t>(globalNode.part)]++;
                    }
                    else {
                        // Node is a ghost.
                        globalNode.globalIdx = nodeGlobalGhostIdxCount++;
                    }
                }  // Finished with this node.
            }
        }  // Finished with all nodes on tile.
    }      // Finished with all tiles.

    ATLAS_ASSERT(nodeGlobalOwnedIdxCount == nNodesUnique + 1);
    ATLAS_ASSERT(nodeGlobalGhostIdxCount == nNodesTotal + 1);
    ATLAS_ASSERT(idxSum(nodeLocalIdxCount) == nNodesUnique);

    // ---------------------------------------------------------------------------
    // 4. LOCATE LOCAL CELLS.
    //    Now that we know where all the nodes and cells are, we can make a list
    //    of local cells to add the the mesh. "Owner" cells correspond to grid
    //    points on this parition. "Halo" cells are at most nHalo grid points
    //    away from an owner cell.
    // ---------------------------------------------------------------------------

    // Make vector of local cells.
    auto localCells = std::vector<LocalElem>{};

    // Loop over all possible local cells.
    for (idx_t t = 0; t < 6; ++t) {
        // Limit range to bounds recorded earlier.
        const BoundingBox& bounds = cellBounds[static_cast<size_t>(t)];

        for (idx_t j = bounds.jBegin; j < bounds.jEnd; ++j) {
            for (idx_t i = bounds.iBegin; i < bounds.iEnd; ++i) {
                if (invalidCell(i, j)) {
                    continue;
                }

                // Get cell
                GlobalElem& globalCell = globalCells[getCellIdx(i, j, t)];

                // Check if cell is an owner.
                if (globalCell.part == static_cast<idx_t>(thisPart)) {
                    // Cell is an owner.
                    localCells.emplace_back();
                    LocalElem& localCell = localCells.back();

                    localCell.type      = ElemType::OWNER;
                    localCell.halo      = 0;
                    localCell.t         = t;
                    localCell.i         = i;
                    localCell.j         = j;
                    localCell.globalPtr = &globalCell;
                }
                else {
                    // Cell is halo if there are nearby owners.
                    bool ownerFound = false;
                    idx_t halo      = nHalo;
                    for (idx_t jHalo = j - nHalo; jHalo < j + nHalo + 1; ++jHalo) {
                        for (idx_t iHalo = i - nHalo; iHalo < i + nHalo + 1; ++iHalo) {
                            if (invalidCell(iHalo, jHalo) || (iHalo == i && jHalo == j)) {
                                continue;
                            }

                            // Is there a nearby owner cell?
                            const GlobalElem& ownerCell = globalCells[getCellIdx(iHalo, jHalo, t)];

                            const bool isOwner = ownerCell.part == static_cast<idx_t>(thisPart);

                            ownerFound = ownerFound || isOwner;

                            if (isOwner) {
                                // Determine halo level from cell distance (l-infinity norm).
                                const idx_t dist = std::max(std::abs(iHalo - i), std::abs(jHalo - j));

                                halo = std::min(dist, halo);
                            }
                        }
                    }

                    if (ownerFound) {
                        // Cell is a halo.
                        localCells.emplace_back();
                        LocalElem& localCell = localCells.back();

                        localCell.type      = ElemType::HALO;
                        localCell.halo      = halo;
                        localCell.t         = t;
                        localCell.i         = i;
                        localCell.j         = j;
                        localCell.globalPtr = &globalCell;

                    }  // Finished with halo search.
                }      // Finished with cell.
            }
        }  // Finished with all cells on tile.
    }      // Finished with all tiles.

    // Partition by cell type.
    auto haloBeginIt =
        std::stable_partition(localCells.begin(), localCells.end(),
                              [](const LocalElem& localCell) -> bool { return localCell.type == ElemType::OWNER; });

    // Sort by halo.
    std::stable_sort(haloBeginIt, localCells.end(),
                     [](const LocalElem& cellA, const LocalElem& cellB) -> bool { return cellA.halo < cellB.halo; });

    // Point global cell to local cell. This is needed to determine node halos.
    // Need to determine remote index and partition if we haven't done so already.
    for (LocalElem& localCell : localCells) {
        localCell.globalPtr->localPtr = &localCell;

        if (localCell.globalPtr->remoteIdx == undefinedIdx) {
            const PointTIJ tijGlobal =
                jacobian.ijLocalToGlobal(PointIJ(localCell.i + 0.5, localCell.j + 0.5), localCell.t);

            const idx_t& t = tijGlobal.t();
            const idx_t& i = tijGlobal.ij().iCell();
            const idx_t& j = tijGlobal.ij().jCell();

            const GlobalElem& ownerCell = globalCells[getCellIdx(i, j, t)];

            localCell.globalPtr->remoteIdx = ownerCell.remoteIdx;
            localCell.globalPtr->part      = ownerCell.part;
        }
    }

    // ---------------------------------------------------------------------------
    // 5. LOCATE LOCAL NODES.
    //    We can now locate the nodes surrounding the owner and halo cells on the
    //    mesh.
    // ---------------------------------------------------------------------------

    // Make vector of local nodes.
    auto localNodes = std::vector<LocalElem>{};

    // Loop over all possible local nodes.
    for (idx_t t = 0; t < 6; ++t) {
        // Limit range to bounds recorded earlier.
        const BoundingBox& bounds = cellBounds[static_cast<size_t>(t)];

        for (idx_t j = bounds.jBegin; j < bounds.jEnd + 1; ++j) {
            for (idx_t i = bounds.iBegin; i < bounds.iEnd + 1; ++i) {
                if (invalidNode(i, j)) {
                    continue;
                }

                // Get node.
                GlobalElem& globalNode = globalNodes[getNodeIdx(i, j, t)];

                // Check if node is an owner.
                if (globalNode.part == static_cast<idx_t>(thisPart)) {
                    // Node is an owner.
                    localNodes.emplace_back();
                    LocalElem& localNode = localNodes.back();

                    localNode.type      = ElemType::OWNER;
                    localNode.halo      = 0;
                    localNode.i         = i;
                    localNode.j         = j;
                    localNode.t         = t;
                    localNode.globalPtr = &globalNode;
                }
                else {
                    // Node is possibly a ghost or halo.

                    // Search four neighbouring cells of node.
                    bool isGhost = false;
                    idx_t halo   = nHalo;

                    // A double-loop with no more than four iterations!
                    for (idx_t iCell = i - 1; iCell < i + 1; ++iCell) {
                        for (idx_t jCell = j - 1; jCell < j + 1; ++jCell) {
                            if (invalidCell(iCell, jCell)) {
                                continue;
                            }

                            const GlobalElem& globalCell = globalCells[getCellIdx(iCell, jCell, t)];

                            // Get halo information from cell.
                            // By this point, any global cell on this PE will have a
                            // non-null local cell pointer.
                            if (globalCell.localPtr) {
                                isGhost = true;
                                halo    = std::min<idx_t>(halo, globalCell.localPtr->halo);
                            }
                        }
                    }

                    if (isGhost) {
                        // Node is an ghost/halo.
                        localNodes.emplace_back();
                        LocalElem& localNode = localNodes.back();

                        localNode.type      = ElemType::HALO;
                        localNode.halo      = halo;
                        localNode.t         = t;
                        localNode.i         = i;
                        localNode.j         = j;
                        localNode.globalPtr = &globalNode;
                    }

                }  // Finished with this node.
            }
        }  // Finished with all nodes on tile.
    }      // Finished with all tiles.

    // Partition by node type.
    auto ghostBeginIt =
        std::stable_partition(localNodes.begin(), localNodes.end(),
                              [](const LocalElem& localNode) -> bool { return localNode.type == ElemType::OWNER; });

    // Sort by halo.
    std::stable_sort(ghostBeginIt, localNodes.end(),
                     [](const LocalElem& nodeA, const LocalElem& nodeB) -> bool { return nodeA.halo < nodeB.halo; });

    // Associate global nodes with local nodes. Allows us to determine node
    // local indices around a cell.
    // Determine partition and remote index if we haven't done so already.
    for (LocalElem& localNode : localNodes) {
        localNode.globalPtr->localPtr = &localNode;

        if (localNode.globalPtr->remoteIdx == undefinedIdx) {
            const PointTIJ tijGlobal = jacobian.ijLocalToGlobal(PointIJ(localNode.i, localNode.j), localNode.t);

            const idx_t& t = tijGlobal.t();
            const idx_t& i = tijGlobal.ij().iNode();
            const idx_t& j = tijGlobal.ij().jNode();

            const GlobalElem& ownerNode = globalNodes[getNodeIdx(i, j, t)];

            localNode.globalPtr->remoteIdx = ownerNode.remoteIdx;
            localNode.globalPtr->part      = ownerNode.part;
        }
    }

    // ---------------------------------------------------------------------------
    // 6. ASSIGN NODES TO MESH
    //    In addition to the usual fields, we will add t, i and j to mesh.nodes().
    // ---------------------------------------------------------------------------

    // Resize nodes.
    auto& nodes = mesh.nodes();
    nodes.resize(static_cast<idx_t>(localNodes.size()));

    // Add extra field.
    Field tijField = nodes.add(Field("tij", make_datatype<idx_t>(), make_shape(nodes.size(), 3)));
    tijField.set_variables(3);

    // Get field views
    auto nodesGlobalIdx = array::make_view<gidx_t, 1>(nodes.global_index());
    auto nodesRemoteIdx = array::make_indexview<idx_t, 1>(nodes.remote_index());
    auto nodesXy        = array::make_view<double, 2>(nodes.xy());
    auto nodesLonLat    = array::make_view<double, 2>(nodes.lonlat());
    auto nodesPart      = array::make_view<int, 1>(nodes.partition());
    auto nodesGhost     = array::make_view<int, 1>(nodes.ghost());
    auto nodesHalo      = array::make_view<int, 1>(nodes.halo());
    auto nodesFlags     = array::make_view<int, 1>(nodes.flags());
    auto nodesTij       = array::make_view<idx_t, 2>(tijField);

    // Set fields.
    idx_t nodeLocalIdx = 0;
    for (const LocalElem& localNode : localNodes) {
        // Set global index.
        nodesGlobalIdx(nodeLocalIdx) = localNode.globalPtr->globalIdx;

        // Set node remote index.
        nodesRemoteIdx(nodeLocalIdx) = localNode.globalPtr->remoteIdx;

        // Set node partition.
        nodesPart(nodeLocalIdx) = localNode.globalPtr->part;

        // Set xy.
        const PointXY xyLocal = jacobian.xy(PointIJ(localNode.i, localNode.j), localNode.t);
        const PointXY xyGlobal =
            interiorNode(localNode.i, localNode.j) ? xyLocal : jacobian.xyLocalToGlobal(xyLocal, localNode.t).xy();

        nodesXy(nodeLocalIdx, XX) = xyLocal.x();
        nodesXy(nodeLocalIdx, YY) = xyLocal.y();

        // Set lon-lat.
        const PointLonLat lonLat       = csProjection.lonlat(xyGlobal);
        nodesLonLat(nodeLocalIdx, LON) = lonLat.lon();
        nodesLonLat(nodeLocalIdx, LAT) = lonLat.lat();

        // Set tij.
        nodesTij(nodeLocalIdx, Coordinates::T) = localNode.t;
        nodesTij(nodeLocalIdx, Coordinates::I) = localNode.i;
        nodesTij(nodeLocalIdx, Coordinates::J) = localNode.j;

        // Set halo.
        nodesHalo(nodeLocalIdx) = localNode.halo;

        // Set flags.
        Topology::reset(nodesFlags(nodeLocalIdx));
        switch (localNode.type) {
            case ElemType::UNDEFINED: {
                ATLAS_ASSERT(0, "Trying to add undefined node to mesh.");
                break;
            }
            case ElemType::OWNER: {
                nodesGhost(nodeLocalIdx) = 0;
                // Vitally important that these two match!
                ATLAS_ASSERT(nodeLocalIdx == localNode.globalPtr->remoteIdx,
                             "Owner local index and remote index do not match.");
                break;
            }
            case ElemType::HALO: {
                nodesGhost(nodeLocalIdx) = 1;
                Topology::set(nodesFlags(nodeLocalIdx), Topology::GHOST);
                break;
            }
        }

        ++nodeLocalIdx;
    }

    // ---------------------------------------------------------------------------
    // 7. ASSIGN CELLS TO MESH
    //    Again, we'll add t, i and j to mesh.cells(). We'll also add the xy and
    //    lonlat coordinates of the cell centres.
    // ---------------------------------------------------------------------------

    auto& cells = mesh.cells();

    // Resize cells.
    cells.add(mesh::ElementType::create("Quadrilateral"), static_cast<idx_t>(localCells.size()));

    // Add extra fields.
    tijField = cells.add(Field("tij", make_datatype<idx_t>(), make_shape(cells.size(), 3)));
    tijField.set_variables(3);

    Field xyField = cells.add(Field("xy", make_datatype<double>(), make_shape(cells.size(), 2)));
    xyField.set_variables(2);

    Field lonLatField = cells.add(Field("lonlat", make_datatype<double>(), make_shape(cells.size(), 2)));
    lonLatField.set_variables(2);

    Field ghostField = cells.add(Field("ghost", make_datatype<int>(), make_shape(cells.size())));

    // Set field views.
    auto cellsGlobalIdx = array::make_view<gidx_t, 1>(cells.global_index());
    auto cellsRemoteIdx = array::make_indexview<idx_t, 1>(cells.remote_index());
    auto cellsPart      = array::make_view<int, 1>(cells.partition());
    auto cellsGhost     = array::make_view<int, 1>(ghostField);
    auto cellsHalo      = array::make_view<int, 1>(cells.halo());
    auto cellsFlags     = array::make_view<int, 1>(cells.flags());
    auto cellsTij       = array::make_view<idx_t, 2>(tijField);
    auto cellsXy        = array::make_view<double, 2>(xyField);
    auto cellsLonLat    = array::make_view<double, 2>(lonLatField);

    // Set local cells.
    auto& nodeConnectivity   = cells.node_connectivity();
    const idx_t cellElemIdx0 = cells.elements(0).begin();

    // Set method to get node local index.
    const auto getNodeLocalIdx = [&](idx_t i, idx_t j, idx_t t) -> idx_t {
        const GlobalElem& globalNode = globalNodes[getNodeIdx(i, j, t)];

        // Do some pointer arithmetic to get local index.
        return static_cast<idx_t>(globalNode.localPtr - localNodes.data());
    };

    idx_t cellLocalIdx = 0;
    for (const LocalElem& localCell : localCells) {
        // Get local indices four surroundings nodes.
        const auto quadNodeIdx = std::array<idx_t, 4>{getNodeLocalIdx(localCell.i, localCell.j, localCell.t),
                                                      getNodeLocalIdx(localCell.i + 1, localCell.j, localCell.t),
                                                      getNodeLocalIdx(localCell.i + 1, localCell.j + 1, localCell.t),
                                                      getNodeLocalIdx(localCell.i, localCell.j + 1, localCell.t)};

        // Set connectivity.
        nodeConnectivity.set(cellLocalIdx + cellElemIdx0, quadNodeIdx.data());

        // Set global index.
        cellsGlobalIdx(cellLocalIdx) = localCell.globalPtr->globalIdx;

        // Set cell remote index.
        cellsRemoteIdx(cellLocalIdx) = localCell.globalPtr->remoteIdx;

        // Set partition.
        cellsPart(cellLocalIdx) = localCell.globalPtr->part;

        // Set xy.
        const PointXY xyLocal = jacobian.xy(PointIJ(localCell.i + 0.5, localCell.j + 0.5), localCell.t);
        const PointXY xyGlobal =
            interiorCell(localCell.i, localCell.j) ? xyLocal : jacobian.xyLocalToGlobal(xyLocal, localCell.t).xy();

        cellsXy(cellLocalIdx, XX) = xyLocal.x();
        cellsXy(cellLocalIdx, YY) = xyLocal.y();

        // Set lon-lat.
        const PointLonLat lonLat       = csProjection.lonlat(xyGlobal);
        cellsLonLat(cellLocalIdx, LON) = lonLat.lon();
        cellsLonLat(cellLocalIdx, LAT) = lonLat.lat();

        // Set tij.
        cellsTij(cellLocalIdx, Coordinates::T) = localCell.t;
        cellsTij(cellLocalIdx, Coordinates::I) = localCell.i;
        cellsTij(cellLocalIdx, Coordinates::J) = localCell.j;


        // Set halo.
        cellsHalo(cellLocalIdx) = localCell.halo;

        // Set ghost.
        cellsGhost(cellLocalIdx) = cellsHalo(cellLocalIdx) > 0;

        // Set flags.
        Topology::reset(cellsFlags(cellLocalIdx));
        switch (localCell.type) {
            case ElemType::UNDEFINED: {
                ATLAS_ASSERT(0, "Trying to add undefined cell to mesh.");
                break;
            }
            case ElemType::OWNER: {
                // Vitally important that these two match!
                ATLAS_ASSERT(cellLocalIdx == localCell.globalPtr->remoteIdx,
                             "Owner local index and remote index do not match.");
                break;
            }
            case ElemType::HALO: {
                Topology::set(cellsFlags(cellLocalIdx), Topology::GHOST);
                break;
            }
        }

        ++cellLocalIdx;
    }

    // ---------------------------------------------------------------------------
    // 8. FINALISE
    //    Done. That was rather a lot of bookkeeping!
    // ---------------------------------------------------------------------------

    set_metadata(mesh);

    return;
}

// -----------------------------------------------------------------------------

void CubedSphereMeshGenerator::set_metadata(Mesh& mesh) const {
    const auto nHalo = options.get<int>("halo");

    // Set basic halo metadata.
    mesh.metadata().set("halo", nHalo);
    mesh.metadata().set("halo_locked", true);
    mesh.nodes().metadata().set("parallel", true);
    mesh.cells().metadata().set("parallel", true);

    // Loop over nodes and count number of halo elements.
    auto nNodes         = std::vector<idx_t>(nHalo + 1, 0);
    const auto nodeHalo = array::make_view<int, 1>(mesh.nodes().halo());
    for (idx_t i = 0; i < mesh.nodes().size(); ++i) {
        ++nNodes[static_cast<size_t>(nodeHalo(i))];
    }
    std::partial_sum(nNodes.begin(), nNodes.end(), nNodes.begin());

    // Set node halo metadata.
    for (size_t i = 0; i < nNodes.size(); ++i) {
        const auto str = "nb_nodes_including_halo[" + std::to_string(i) + "]";
        mesh.metadata().set(str, nNodes[i]);
    }
    // Loop over cells and count number of halo elements.
    auto nCells         = std::vector<idx_t>(nHalo + 1, 0);
    const auto cellHalo = array::make_view<int, 1>(mesh.cells().halo());
    for (idx_t i = 0; i < mesh.cells().size(); ++i) {
        ++nCells[static_cast<size_t>(cellHalo(i))];
    }
    std::partial_sum(nCells.begin(), nCells.end(), nCells.begin());

    // Set cell halo metadata.
    for (size_t i = 0; i < nCells.size(); ++i) {
        const auto str = "nb_cells_including_halo[0][" + std::to_string(i) + "]";
        mesh.metadata().set(str, nCells[i]);
    }
}

// -----------------------------------------------------------------------------

void CubedSphereMeshGenerator::hash(eckit::Hash& h) const {
    h.add("CubedSphereMeshGenerator");
    options.hash(h);
}

// -----------------------------------------------------------------------------

namespace {
static MeshGeneratorBuilder<CubedSphereMeshGenerator> CubedSphereMeshGenerator(CubedSphereMeshGenerator::static_type());
}

// -----------------------------------------------------------------------------

}  // namespace meshgenerator
}  // namespace atlas
