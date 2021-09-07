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

// -----------------------------------------------------------------------------

CubedSphereMeshGenerator::CubedSphereMeshGenerator( const eckit::Parametrisation& p ) {
    configure_defaults();

    // Get number of partitions.
    size_t nb_parts;
    if ( p.get( "nb_parts", nb_parts ) )
        options.set( "nb_parts", nb_parts );

    // Get this partition.
    size_t part;
    if ( p.get( "part", part ) )
        options.set( "part", part );

    // Get partitioner.
    std::string partitioner;
    if ( p.get( "partitioner", partitioner ) )
        options.set( "partitioner", partitioner );
}

// -----------------------------------------------------------------------------


void CubedSphereMeshGenerator::configure_defaults() {
    // This option sets number of parts the mesh will be split in.
    options.set( "nb_parts", mpi::size() );

    // This option sets the part that will be generated.
    options.set( "part", mpi::rank() );

    // This options sets the default partitioner.
    options.set<std::string>( "partitioner", "cubedsphere" );
}

// -----------------------------------------------------------------------------

void CubedSphereMeshGenerator::generate( const Grid& grid, Mesh& mesh ) const {
    // Get partitioner type and number of partitions from config.
    const auto nParts   = static_cast<idx_t>( options.get<size_t>( "nb_parts" ) );
    const auto partType = options.get<std::string>( "partitioner" );
    auto partConfig     = util::Config{};
    partConfig.set( "type", partType );
    partConfig.set( "partitions", nParts );

    // Use lonlat instead of xy for equal_regions partitioner.
    if ( partType == "equal_regions" )
        partConfig.set( "coordinates", "lonlat" );

    // Set distribution.
    const auto partitioner  = grid::Partitioner( partConfig );
    const auto distribution = grid::Distribution( grid, partitioner );

    generate( grid, distribution, mesh );
}

// -----------------------------------------------------------------------------

void CubedSphereMeshGenerator::generate( const Grid& grid, const grid::Distribution& distribution, Mesh& mesh ) const {
    // Check for correct grid and need for mesh
    ATLAS_ASSERT( !mesh.generated() );
    if ( !CubedSphereGrid( grid ) ) {
        throw_Exception(
            "CubedSphereMeshGenerator can only work "
            "with a cubedsphere grid.",
            Here() );
    }

    // Cast grid to cubed sphere grid.
    const auto csGrid = CubedSphereGrid( grid );

    // Check for correct stagger
    const auto gridName = grid.name();

    if ( csGrid.stagger() != "C" ) {
        throw_Exception(
            "CubedSphereMeshGenerator can only work with a "
            "cell-centroid grid. Try NodalCubedSphereMeshGenerator instead." );
    }


    // Clone some grid properties.
    setGrid( mesh, csGrid, distribution );

    generate_mesh( csGrid, distribution, mesh );
}

// -----------------------------------------------------------------------------

void CubedSphereMeshGenerator::generate_mesh( const CubedSphereGrid& csGrid, const grid::Distribution& distribution,
                                              Mesh& mesh ) const {
    ATLAS_TRACE( "CubedSphereMeshGenerator::generate" );

    // ---------------------------------------------------------------------------
    // CUBED SPHERE MESH GENERATOR
    // ---------------------------------------------------------------------------
    //
    // Mesh generator creates a mesh by placing cells with centroid positions at
    // grid xy coordinates. Nodes are then placed around cells by the generator.
    //
    // The node-ownership of cells is determined by the following rules:
    //    * node (i > 0, j > 0) is owned by cell (i - 1, j - 1)
    //    * node (i = 0, j > 0) is owned by cell (i    , j - 1)
    //    * node (i > 0, j = 0) is owned by cell (i - 1, j    )
    //    * node (i = 0, j = 0) is owned by cell (i    , j    )
    //
    // In addtion, there are two types of ghost node that need to be added to the
    // mesh:
    //    * Type-A ghost nodes. These are added to an edge between two tiles. One
    //      tile owns the nodes, the other tile has ghost nodes in the same place.
    //      The ownership rules are determined the CubedSphereTiles class. These
    //      nodes are defined globally and have their own unique global index.
    //      Their placement and local indexing is unaffected by the partitioner.
    //    * Type-B ghost nodes. These are added when a tile is partitioned. One
    //      side of the partation boundary has cells which own the nodes
    //      (see ownership rules above). The other side has ghost nodes. The
    //      placement and indexing of these ghost nodes may vary with different
    //      partitioners. The global index of these nodes are taken from their
    //      owned counterparts on the other side of the partition boundary.
    //
    // Global indices of cells are the same as the grid point indices. Global
    // indices of nodes first count the owned nodes, then the type-A ghost points.
    //
    // Local indices of cells follows the order of the local grid points. Local
    // indces of the nodes first count the owned nodes, followed by type-A ghost
    // nodes, followed by type-B ghost nodes.
    //
    // There are several stages to the mesh generator:
    //    1. Preamble.
    //    2. Define global cell distribution.
    //    3. Define global owned node and type-A ghost node distribution.
    //    4. Locate and count local type-B ghost nodes.
    //    5. Assign nodes to mesh.
    //    6. Assign cells to mesh.
    //    7. Finalise.

    // ---------------------------------------------------------------------------
    // 1. PREAMBLE
    // ---------------------------------------------------------------------------

    // Get dimensions of grid
    const idx_t N      = csGrid.N();
    const idx_t nTiles = csGrid.tiles().size();

    const auto nNodesUnique = nTiles * N * N + 2;
    const auto nNodesAll    = nTiles * ( N + 1 ) * ( N + 1 );
    const auto nCells       = nTiles * N * N;

    Log::debug() << "Number of cells per tile edge = " << std::to_string( N ) << std::endl;

    // Define bad index values.
    constexpr auto undefinedIdx       = -1;
    constexpr auto undefinedGlobalIdx = -1;

    // Projection and tiles
    const auto* const csProjection =
        dynamic_cast<const projection::detail::CubedSphereProjectionBase*>( csGrid.projection().get() );
    ATLAS_ASSERT( csProjection );
    const auto csTiles = csProjection->getCubedSphereTiles();

    // Get partition information.
    const auto nParts   = mpi::comm().size();
    const auto thisPart = mpi::comm().rank();

    // Set some casting helper functions to avoid annoying warnings.
    const auto st2idx = []( size_t i ) { return static_cast<idx_t>( i ); };
    const auto idx2st = []( idx_t i ) { return static_cast<size_t>( i ); };

    // Define an index counter.
    const auto idxSum = []( const std::vector<idx_t>& idxCounts ) -> idx_t {
        return std::accumulate( idxCounts.begin(), idxCounts.end(), 0 );
    };

    // Helper functions to get node and cell idx from (t, j, i).
    const auto getNodeIdx = [&]( idx_t t, idx_t j, idx_t i ) {
        // Adjust bounds.
        constexpr idx_t zero = 0;
        t                    = std::max( std::min<idx_t>( t, nTiles - 1 ), zero );
        j                    = std::max( std::min( j, N ), zero );
        i                    = std::max( std::min( i, N ), zero );
        return idx2st( t * ( N + 1 ) * ( N + 1 ) + j * ( N + 1 ) + i );
    };

    const auto getCellIdx = [&]( idx_t t, idx_t j, idx_t i ) {
        // Adjust bounds.
        constexpr idx_t zero = 0;
        t                    = std::max( std::min( t, nTiles - 1 ), zero );
        j                    = std::max( std::min( j, N - 1 ), zero );
        i                    = std::max( std::min( i, N - 1 ), zero );
        return idx2st( t * N * N + j * N + i );
    };

    // ---------------------------------------------------------------------------
    // 2. GLOBAL CELL DISTRIBUTION
    // ---------------------------------------------------------------------------

    // Define cell record.
    // This currently keeps track of more information than we probably need.
    // This may be reduced once we've implemented halos.
    struct CellRecord {
        gidx_t globalIdx{undefinedGlobalIdx};  // Global ID.
        idx_t localIdx{undefinedIdx};          // Local ID.
        idx_t part{undefinedIdx};              // Partition.
        PointXY xy{};                          // Position.
    };

    // Define ij bounding box for each face (this partition only).
    struct BoundingBox {
        idx_t iBegin{std::numeric_limits<idx_t>::max()};
        idx_t iEnd{std::numeric_limits<idx_t>::min()};
        idx_t jBegin{std::numeric_limits<idx_t>::max()};
        idx_t jEnd{std::numeric_limits<idx_t>::min()};
    };

    // Make list of all cells.
    auto globalCells = std::vector<CellRecord>( idx2st( nCells ) );

    // Initialise bounding box.
    auto cellBounds = std::vector<BoundingBox>( idx2st( nTiles ) );

    // Set xy and tji grid iterators.
    auto tjiIt = csGrid.tij().begin();
    auto xyIt  = csGrid.xy().begin();

    // Set counters for cell local indices.
    auto cellIdxCount = std::vector<idx_t>( nParts, 0 );

    for ( gidx_t gridIdx = 0; gridIdx < csGrid.size(); ++gridIdx ) {
        // Grid (t, j, i) order does not have to match mesh (t, j, i), although
        // in practice it probably does.

        // Get cell index.
        const auto t       = ( *tjiIt ).t();
        const auto j       = ( *tjiIt ).j();
        const auto i       = ( *tjiIt ).i();
        const auto cellIdx = getCellIdx( t, j, i );
        auto& cell         = globalCells[cellIdx];

        // Set global index.
        cell.globalIdx = gridIdx + 1;

        // Set partition and remote index.
        cell.part = distribution.partition( gridIdx );

        // Set local index.
        cell.localIdx = cellIdxCount[idx2st( cell.part )]++;

        // Set xy
        cell.xy = *xyIt;

        if ( cell.part == st2idx( thisPart ) ) {
            // Keep track of local (t, j, i) bounds.
            auto& bounds  = cellBounds[idx2st( t )];
            bounds.iBegin = std::min( bounds.iBegin, i );
            bounds.iEnd   = std::max( bounds.iEnd, i + 1 );
            bounds.jBegin = std::min( bounds.jBegin, j );
            bounds.jEnd   = std::max( bounds.jEnd, j + 1 );
        }
        // Increment iterators.
        ++tjiIt;
        ++xyIt;
    }

    ATLAS_ASSERT( idxSum( cellIdxCount ) == nCells );

    // ---------------------------------------------------------------------------
    // 3. GLOBAL OWNED NODE AND TYPE-A GHOST NODE DISTRIBUTION
    // ---------------------------------------------------------------------------

    // Define Jacobian of xy wrt ij for each tile so that we can easily
    // switch between the two. This will be needed to set node xy positions
    // and figure out the ij indices for type-A ghost nodes.
    struct XyJacobian {
        double dxByDi{};  // dx/di
        double dxByDj{};  // dx/dj
        double dyByDi{};  // dy/di
        double dyByDj{};  // dy/dj
        double diByDx{};  // di/dx
        double diByDy{};  // di/dy
        double djByDx{};  // dj/dx
        double djByDy{};  // dj/dy
        PointXY xy00{};   // xy coordinate of node(i = 0, j = 0).
    };

    // Calculate Jacobians.
    idx_t t   = 0;
    auto jacs = std::vector<XyJacobian>{};
    std::generate_n( std::back_inserter( jacs ), nTiles, [&]() {
        // Initialise element.
        auto jacElem = XyJacobian{};

        // Get cell positions.
        const auto& xy00 = globalCells[getCellIdx( t, 0, 0 )].xy;
        const auto& xy10 = globalCells[getCellIdx( t, 0, 1 )].xy;
        const auto& xy01 = globalCells[getCellIdx( t, 1, 0 )].xy;

        // Calculate Jacobian.
        jacElem.dxByDi = xy10.x() - xy00.x();
        jacElem.dxByDj = xy01.x() - xy00.x();
        jacElem.dyByDi = xy10.y() - xy00.y();
        jacElem.dyByDj = xy01.y() - xy00.y();

        // Calculate inverse Jacobian.
        const auto invDet = 1. / ( jacElem.dxByDi * jacElem.dyByDj - jacElem.dxByDj * jacElem.dyByDi );
        jacElem.diByDx    = jacElem.dyByDj * invDet;
        jacElem.diByDy    = -jacElem.dxByDj * invDet;
        jacElem.djByDx    = -jacElem.dyByDi * invDet;
        jacElem.djByDy    = jacElem.dxByDi * invDet;

        // Extrapolate cell(t, 0, 0) xy to get node(t, 0, 0) xy.
        jacElem.xy00 = PointXY{xy00.x() - 0.5 * jacElem.dxByDi - 0.5 * jacElem.dxByDj,
                               xy00.y() - 0.5 * jacElem.dyByDi - 0.5 * jacElem.dyByDj};

        ++t;
        return jacElem;
    } );

    // Define node record.
    // This currently keeps track of more information than we probably need.
    // This may be reduced once we've implemented halos.
    struct NodeRecord {
        gidx_t globalIdx{undefinedGlobalIdx};  // Global ID,
        idx_t localIdx{undefinedIdx};          // Local ID.
        idx_t localPart{undefinedIdx};         // Partition node is on.
        idx_t remoteIdx{undefinedIdx};         // Local ID of owned node.
        idx_t remotePart{undefinedIdx};        // Partion that owned node is on.
        idx_t t{undefinedIdx};                 // Tile that owned node is on.
        bool ghost{false};                     // True if node is a ghost.
        PointXY localXy{};                     // Position of this node.
        PointXY remoteXy{};                    // Position of owned node.
    };

    // Make list of all nodes.
    auto globalNodes = std::vector<NodeRecord>( idx2st( nNodesAll ) );

    // Set counters for local node indices.
    auto nodeLocalIdxCount = std::vector<idx_t>( nParts, 0 );

    // Set counter for global indices.
    gidx_t nodeGlobalIdxCount = 1;

    // Make a list of type-A ghost nodes and process later.
    auto typeAGhostNodes = std::vector<NodeRecord*>{};

    for ( idx_t t = 0; t < nTiles; ++t ) {
        for ( idx_t j = 0; j < N + 1; ++j ) {
            for ( idx_t i = 0; i < N + 1; ++i ) {
                // Get node.
                auto& node = globalNodes[getNodeIdx( t, j, i )];

                // Get owning cell.
                const auto& cell = globalCells[getCellIdx( t, j - 1, i - 1 )];

                // Extrapolate local xy from cell centre.
                // Nodes are +/- half a grid spacing relative to cells.
                const auto di   = i == 0 ? -0.5 : 0.5;
                const auto dj   = j == 0 ? -0.5 : 0.5;
                const auto& jac = jacs[idx2st( t )];
                node.localXy    = PointXY{cell.xy.x() + di * jac.dxByDi + dj * jac.dxByDj,
                                       cell.xy.y() + di * jac.dyByDi + dj * jac.dyByDj};

                // Get remote xy from tile class and local t.
                node.remoteXy = csTiles.tileCubePeriodicity( node.localXy, t );

                // Local taken from owning cell.
                node.localPart = cell.part;

                // Get remote t from tile class.
                node.t = csTiles.indexFromXY( node.remoteXy.data() );

                // Node is a type-A ghost if local t and remote t differ.
                // Otherwise node is owned.
                if ( t == node.t ) {
                    // Owned point.

                    // Set global index.
                    node.globalIdx = nodeGlobalIdxCount++;

                    // Set local index.
                    node.localIdx = nodeLocalIdxCount[idx2st( node.localPart )]++;

                    // Remote partition same as local.
                    node.remotePart = node.localPart;

                    // Remote index same as local.
                    node.remoteIdx = node.localIdx;
                }
                else {
                    // Type-A ghost. Deal with this later.
                    typeAGhostNodes.push_back( &node );
                }
            }
        }
    }

    ATLAS_ASSERT( nodeGlobalIdxCount == nNodesUnique + 1 );
    ATLAS_ASSERT( idxSum( nodeLocalIdxCount ) == nNodesUnique );

    // Deal with type-A ghost nodes.
    for ( auto* nodePtr : typeAGhostNodes ) {
        // Set global index.
        nodePtr->globalIdx = nodeGlobalIdxCount++;

        // Set local index.
        nodePtr->localIdx = nodeLocalIdxCount[idx2st( nodePtr->localPart )]++;

        // Get jacobian.
        const auto& jac = jacs[idx2st( nodePtr->t )];

        // Get actual i and j.
        const auto dx = nodePtr->remoteXy.x() - jac.xy00.x();
        const auto dy = nodePtr->remoteXy.y() - jac.xy00.y();
        const auto i  = static_cast<idx_t>( std::round( dx * jac.diByDx + dy * jac.diByDy ) );
        const auto j  = static_cast<idx_t>( std::round( dx * jac.djByDx + dy * jac.djByDy ) );

        // Get remote node.
        const auto& remoteNode = globalNodes[getNodeIdx( nodePtr->t, j, i )];

        // Set partition and remote index.
        nodePtr->remotePart = remoteNode.localPart;
        nodePtr->remoteIdx  = remoteNode.localIdx;

        // Node is a ghost.
        nodePtr->ghost = true;

        // Check that remote index is now defined.
        ATLAS_ASSERT( nodePtr->remoteIdx != undefinedIdx );
    }

    ATLAS_ASSERT( nodeGlobalIdxCount == nNodesAll + 1 );
    ATLAS_ASSERT( idxSum( nodeLocalIdxCount ) == nNodesAll );

    // ---------------------------------------------------------------------------
    // 4. LOCATE AND COUNT LOCAL TYPE-B GHOST NODES.
    // ---------------------------------------------------------------------------

    // Make list of local nodes (this simplifies part 5).
    auto localNodes = std::vector<const NodeRecord*>{};

    // Loop over all possible local nodes.
    for ( idx_t t = 0; t < nTiles; ++t ) {
        // Limit range to bounds recorded earlier.
        const auto bounds = cellBounds[idx2st( t )];

        for ( idx_t j = bounds.jBegin; j < bounds.jEnd + 1; ++j ) {
            for ( idx_t i = bounds.iBegin; i < bounds.iEnd + 1; ++i ) {
                // Get node.
                auto& node = globalNodes[getNodeIdx( t, j, i )];

                // Get four neighbouring cells of node. (cell0 is owner of node).
                const auto& cell0 = globalCells[getCellIdx( t, j - 1, i - 1 )];
                const auto& cell1 = globalCells[getCellIdx( t, j - 1, i )];
                const auto& cell2 = globalCells[getCellIdx( t, j, i )];
                const auto& cell3 = globalCells[getCellIdx( t, j, i - 1 )];

                if ( cell0.part == st2idx( thisPart ) ) {
                    // Check that node-cell partitions match.
                    ATLAS_ASSERT( node.localPart == cell0.part );

                    // Node is either owned or Type-A ghost.
                    localNodes.push_back( &node );
                }
                else if ( cell1.part == st2idx( thisPart ) || cell2.part == st2idx( thisPart ) ||
                          cell3.part == st2idx( thisPart ) ) {
                    // Node has a non-owning neighbouring cell and is a type-B ghost.

                    // Change local index to match this partition.
                    node.localPart = st2idx( thisPart );
                    node.localIdx  = nodeLocalIdxCount[thisPart]++;
                    node.ghost     = true;

                    // Check remote index is defined.
                    ATLAS_ASSERT( node.remoteIdx != undefinedIdx );

                    localNodes.push_back( &node );
                }
            }
        }
    }

    // ---------------------------------------------------------------------------
    // 5. ASSIGN NODES TO MESH
    // ---------------------------------------------------------------------------

    // Resize nodes.
    mesh.nodes().resize( nodeLocalIdxCount[thisPart] );

    // Get field views
    auto nodesGlobalIdx = array::make_view<gidx_t, 1>( mesh.nodes().global_index() );
    auto nodesRemoteIdx = array::make_indexview<idx_t, 1>( mesh.nodes().remote_index() );
    auto nodesXy        = array::make_view<double, 2>( mesh.nodes().xy() );
    auto nodesLonLat    = array::make_view<double, 2>( mesh.nodes().lonlat() );
    auto nodesPart      = array::make_view<int, 1>( mesh.nodes().partition() );
    auto nodesGhost     = array::make_view<int, 1>( mesh.nodes().ghost() );
    auto nodesFlags     = array::make_view<int, 1>( mesh.nodes().flags() );

    // Set fields.
    for ( const auto* nodePtr : localNodes ) {
        // Get local index.
        const auto localIdx = nodePtr->localIdx;

        // Set global index.
        nodesGlobalIdx( localIdx ) = nodePtr->globalIdx;

        // Set node remote index.
        nodesRemoteIdx( localIdx ) = nodePtr->remoteIdx;

        // Set node partition.
        nodesPart( localIdx ) = nodePtr->remotePart;

        // Set xy.
        nodesXy( localIdx, XX ) = nodePtr->localXy.x();
        nodesXy( localIdx, YY ) = nodePtr->localXy.y();

        // Set lon-lat.
        const auto lonLat            = csProjection->lonlat( nodePtr->remoteXy );
        nodesLonLat( localIdx, LON ) = lonLat.lon();
        nodesLonLat( localIdx, LAT ) = lonLat.lat();

        // Set ghost flag.
        nodesGhost( localIdx ) = nodePtr->ghost;

        // Set general flags.
        Topology::reset( nodesFlags( localIdx ) );
        if ( nodePtr->ghost )
            Topology::set( nodesFlags( localIdx ), Topology::GHOST );
    }

    // ---------------------------------------------------------------------------
    // 6. ASSIGN CELLS TO MESH
    // ---------------------------------------------------------------------------

    // Resize cells.
    mesh.cells().add( new mesh::temporary::Quadrilateral(), cellIdxCount[thisPart] );

    // Set field views.
    auto cellsGlobalIdx = array::make_view<gidx_t, 1>( mesh.cells().global_index() );
    auto cellsPart      = array::make_view<int, 1>( mesh.cells().partition() );

    // Set local cells.
    auto& nodeConnectivity = mesh.cells().node_connectivity();
    const auto idx0        = mesh.cells().elements( 0 ).begin();

    for ( idx_t t = 0; t < nTiles; ++t ) {
        // Use bounds from before.
        const auto bounds = cellBounds[idx2st( t )];

        for ( idx_t j = bounds.jBegin; j < bounds.jEnd; ++j ) {
            for ( idx_t i = bounds.iBegin; i < bounds.iEnd; ++i ) {
                // Get cell.
                const auto& cell = globalCells[getCellIdx( t, j, i )];

                // Only add cells on this partition.
                if ( cell.part == st2idx( thisPart ) ) {
                    // Get local index.
                    const auto localIdx = cell.localIdx + idx0;

                    // Set quadrilateral.
                    const auto& node0 = globalNodes[getNodeIdx( t, j, i )];
                    const auto& node1 = globalNodes[getNodeIdx( t, j, i + 1 )];
                    const auto& node2 = globalNodes[getNodeIdx( t, j + 1, i + 1 )];
                    const auto& node3 = globalNodes[getNodeIdx( t, j + 1, i )];

                    // Check nodes partitions match.
                    ATLAS_ASSERT( node0.localPart == cell.part );
                    ATLAS_ASSERT( node1.localPart == cell.part );
                    ATLAS_ASSERT( node2.localPart == cell.part );
                    ATLAS_ASSERT( node3.localPart == cell.part );

                    const auto quadNodeIdx =
                        std::array<idx_t, 4>{node0.localIdx, node1.localIdx, node2.localIdx, node3.localIdx};

                    // Set connectivity.
                    nodeConnectivity.set( localIdx, quadNodeIdx.data() );

                    // Set global index.
                    cellsGlobalIdx( localIdx ) = cell.globalIdx;

                    // Set partition.
                    cellsPart( localIdx ) = cell.part;
                }
            }
        }
    }

    // ---------------------------------------------------------------------------
    // 7. FINALISE
    // ---------------------------------------------------------------------------

    mesh.nodes().global_index().metadata().set( "human_readable", true );
    mesh.nodes().global_index().metadata().set( "min", 1 );
    mesh.nodes().global_index().metadata().set( "max", nNodesAll );
    mesh.nodes().metadata().set( "parallel", true );

    mesh.cells().global_index().metadata().set( "human_readable", true );
    mesh.cells().global_index().metadata().set( "min", 1 );
    mesh.cells().global_index().metadata().set( "max", nCells );
    mesh.cells().metadata().set( "parallel", true );

    return;
}

// -----------------------------------------------------------------------------

void CubedSphereMeshGenerator::hash( eckit::Hash& h ) const {
    h.add( "CubedSphereMeshGenerator" );
    options.hash( h );
}

// -----------------------------------------------------------------------------

namespace {
static MeshGeneratorBuilder<CubedSphereMeshGenerator> CubedSphereMeshGenerator(
    CubedSphereMeshGenerator::static_type() );
}

// -----------------------------------------------------------------------------

}  // namespace meshgenerator
}  // namespace atlas
