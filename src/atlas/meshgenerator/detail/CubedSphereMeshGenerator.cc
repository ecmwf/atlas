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

CubedSphereMeshGenerator::CubedSphereMeshGenerator( const eckit::Parametrisation& p ) {}

// -------------------------------------------------------------------------------------------------

void CubedSphereMeshGenerator::configure_defaults() {}

// -------------------------------------------------------------------------------------------------

void CubedSphereMeshGenerator::generate( const Grid& grid, Mesh& mesh ) const {
    // Check for proper grid and need for mesh
    ATLAS_ASSERT( !mesh.generated() );
    const CubedSphereGrid csg = CubedSphereGrid( grid );
    if ( !csg ) {
        throw_Exception( "CubedSphereMeshGenerator can only work with a cubedsphere grid", Here() );
    }

    // Number of processors
    //size_t nb_parts = options.get<size_t>( "nb_parts" );

    // Decomposition type
    std::string partitioner_type = "checkerboard";
    options.get( "checkerboard", partitioner_type );

    // Partitioner
    grid::Partitioner partitioner( partitioner_type, 1 );
    grid::Distribution distribution( partitioner.partition( grid ) );

    generate( grid, distribution, mesh );
}

// -------------------------------------------------------------------------------------------------

void CubedSphereMeshGenerator::generate( const Grid& grid, const grid::Distribution& distribution, Mesh& mesh ) const {
    const auto csgrid = CubedSphereGrid( grid );

    const int N      = csgrid.N();
    const int nTiles = csgrid.GetNTiles();


    // Make a list linking ghost (t, i, j) values to known (t, i, j)
    // This will be different for each cube-sphere once MeshGenerator is generalised.

    using Tij = typename std::array<idx_t, 3>;

    // first: ghost (t, i, j); second owned (t, i, j).
    using TijPair = typename std::pair<Tij, Tij>;

    auto ghostToOwnedTij = std::vector<TijPair>{};


    // -------------------------------------------------------------------------
    // BEGIN FV3 SPECIFIC MAP
    // -------------------------------------------------------------------------

    // Special points.
    ghostToOwnedTij.push_back(TijPair{{0, N, N}, {2, 0, 0}});
    ghostToOwnedTij.push_back(TijPair{{1, N, N}, {3, 0, 0}});
    ghostToOwnedTij.push_back(TijPair{{2, N, N}, {4, 0, 0}});
    ghostToOwnedTij.push_back(TijPair{{3, N, N}, {5, 0, 0}});
    ghostToOwnedTij.push_back(TijPair{{4, N, N}, {0, 0, 0}});
    ghostToOwnedTij.push_back(TijPair{{5, N, N}, {1, 0, 0}});
    ghostToOwnedTij.push_back(TijPair{{2, 0, N}, {0, 0, N}});
    ghostToOwnedTij.push_back(TijPair{{4, 0, N}, {0, 0, N}});
    ghostToOwnedTij.push_back(TijPair{{3, N, 0}, {1, N, 0}});
    ghostToOwnedTij.push_back(TijPair{{5, N, 0}, {1, N, 0}});

    // Tile 1
    for ( idx_t ix = 1; ix < N; ix++ ) ghostToOwnedTij.push_back(TijPair{{0, ix, N}, {2, 0, N - ix}});
    for ( idx_t iy = 0; iy < N; iy++ ) ghostToOwnedTij.push_back(TijPair{{0, N, iy}, {1, 0, iy}});

    // Tile 2
    for ( idx_t ix = 0; ix < N; ix++ ) ghostToOwnedTij.push_back(TijPair{{1, ix, N}, {2, ix, 0}});
    for ( idx_t iy = 1; iy < N; iy++ ) ghostToOwnedTij.push_back(TijPair{{1, N, iy}, {3, N - iy, 0}});

    // Tile 3
    for ( idx_t ix = 1; ix < N; ix++ ) ghostToOwnedTij.push_back(TijPair{{2, ix, N}, {4, 0, N - ix}});
    for ( idx_t iy = 0; iy < N; iy++ ) ghostToOwnedTij.push_back(TijPair{{2, N, iy}, {3, 0, iy}});

    // Tile 4
    for ( idx_t ix = 0; ix < N; ix++ ) ghostToOwnedTij.push_back(TijPair{{3, ix, N}, {4, ix, 0}});
    for ( idx_t iy = 1; iy < N; iy++ ) ghostToOwnedTij.push_back(TijPair{{3, N, iy}, {5, N - iy, 0}});

    // Tile 5
    for ( idx_t ix = 1; ix < N; ix++ ) ghostToOwnedTij.push_back(TijPair{{4, ix, N}, {0, 0, N - ix}});
    for ( idx_t iy = 0; iy < N; iy++ ) ghostToOwnedTij.push_back(TijPair{{4, N, iy}, {5, 0, iy}});

    // Tile 6
    for ( idx_t ix = 0; ix < N; ix++ ) ghostToOwnedTij.push_back(TijPair{{5, ix, N}, {0, ix, 0}});
    for ( idx_t iy = 1; iy < N; iy++ ) ghostToOwnedTij.push_back(TijPair{{5, N, iy}, {1, N - iy, 0}});



    // -------------------------------------------------------------------------
    // END FV3 SPECIFIC MAP
    // -------------------------------------------------------------------------
    std::cout << " after fv3"  << std::endl;

    ATLAS_TRACE( "CubedSphereMeshGenerator::generate" );
    Log::debug() << "Number of faces per tile edge = " << std::to_string( N ) << std::endl;

    // Number of nodes
    int nnodes     = nTiles * N * N + 2;              // Number of unique grid nodes
    int nnodes_all = nTiles * ( N + 1 ) * ( N + 1 );  // Number of grid nodes including edge and corner duplicates
    int ncells     = nTiles * N * N;                  // Number of unique grid cells

    // Construct mesh nodes
    // --------------------
    mesh.nodes().resize( nnodes_all );
    mesh::Nodes& nodes = mesh.nodes();
    auto xy            = array::make_view<double, 2>( nodes.xy() );
    auto lonlat        = array::make_view<double, 2>( nodes.lonlat() );
    auto glb_idx       = array::make_view<gidx_t, 1>( nodes.global_index() );
    auto remote_idx    = array::make_indexview<idx_t, 1>( nodes.remote_index() );
    auto part          = array::make_view<int, 1>( nodes.partition() );
    auto ghost         = array::make_view<int, 1>( nodes.ghost() );
    auto flags         = array::make_view<int, 1>( nodes.flags() );

    // Loop over entire grid
    // ---------------------

    // Node array will include duplicates of the grid points to make it easier to fill up the
    // neighbours. However the mesh will not contain these points so no Ghost nodes are required.
    // We could include the duplicate points in the array if we make them Ghost nodes but it is not
    // clear if this provides any benefit.



    array::ArrayT<int> NodeArrayT( nTiles, N + 1, N + 1 );  // All grid points including duplicates
    auto NodeArray = array::make_view<int, 3>( NodeArrayT );

    // Add owned nodes to node array
    idx_t nOwned = 0;
    auto addOwnedNode = [&](Tij tijOwned) {

        // Set node array
        NodeArray(tijOwned[0], tijOwned[1], tijOwned[2]) = nOwned;

        // Get xy from global xy grid array
        double xy_[3];
        csgrid.xy( tijOwned[1], tijOwned[2], tijOwned[0], xy_ );

        xy( nOwned, XX ) = xy_[XX];
        xy( nOwned, YY ) = xy_[YY];

        // Get lonlat from global lonlat array
        double lonlat_[2];
        csgrid.lonlat( tijOwned[1], tijOwned[2], tijOwned[0], lonlat_ );

        lonlat( nOwned, LON ) = lonlat_[LON];
        lonlat( nOwned, LAT ) = lonlat_[LAT];

        // Is not ghost node
        mesh::Nodes::Topology::reset(flags(nOwned));
        ghost(nOwned) = 0;

        glb_idx(nOwned) = nOwned + 1;
        remote_idx(nOwned) = nOwned;
        part(nOwned) = distribution.partition(nOwned);

        ++nOwned;

        return;
    };



    auto p  = csgrid.tij().begin();
    auto e = csgrid.tij().end();
    std::cout << "begin  = " << (*p).t()  << " " << (*p).i() << " " << (*p).j() << std::endl;
    std::cout << "end  = " << (*e).t()  << " " << (*e).i() << " " << (*e).j() << std::endl;

    // Loop over owned (t, i, j)
   // for (auto& p : csgrid.tij()) {
    for ( std::size_t jn = 0; jn < csgrid.size(); ++jn ) {
        std::cout << "add Owned Node = " << jn << " " << (*p).t()  << " " << (*p).i() << " " << (*p).j() << std::endl;
        addOwnedNode(Tij{(*p).t(), (*p).i(), (*p).j()});
        if (jn != static_cast<std::size_t>(csgrid.size() - 1)) {++p;}

   }

    // Assert that the correct number of nodes have been set
    ATLAS_ASSERT( nnodes == nOwned, "Insufficient nodes" );

    // Vector of ghost global index of each ghost point
    auto ghostGblIdx = std::vector<idx_t>();
    auto ownedGblIdx = std::vector<idx_t>();

    std::cout << "ghostToOwnedTij.size() " << ghostToOwnedTij.size() << std::endl;


    // Add ghost nodes to node array
    // (nGhost started after nOwned)
    auto nGhost = nOwned;
    auto addGhostNode = [&](Tij tijGhost, Tij tijOwned) {

      std::cout << "addGhostNode" << std::endl;

      // Get concrete node id
      auto nOwned = NodeArray(tijOwned[0], tijOwned[1], tijOwned[2]);

      std::cout << "addGhostNode nOwned id" << nOwned << std::endl;

      // Add ghost node to NodeArray
      NodeArray(tijGhost[0], tijGhost[1], tijGhost[2]) = nGhost;

     std::cout << "addGhostNode to NodeArray" << nOwned << std::endl;

      // "Create" ghost xy coordinate.

      // Get Jacobian of coords rtw indices
      auto x0 = xy(NodeArray(tijGhost[0], 0, 0), XX);

      std::cout << "addGhostNode x0" << x0 << std::endl;

      auto y0 = xy(NodeArray(tijGhost[0], 0, 0), YY);

      std::cout << "addGhostNode y0" << y0 << std::endl;

      auto dx_di = xy(NodeArray(tijGhost[0], 1, 0), XX) - x0;
      auto dx_dj = xy(NodeArray(tijGhost[0], 0, 1), XX) - x0;
      auto dy_di = xy(NodeArray(tijGhost[0], 1, 0), YY) - y0;
      auto dy_dj = xy(NodeArray(tijGhost[0], 0, 1), YY) - y0;

      std::cout << "addGhostNode dx_ " << dx_di << " " <<dx_dj << std::endl;

      // Set xy coordinates
      xy(nGhost, XX) = x0 + tijGhost[1] * dx_di + tijGhost[2] * dx_dj;
      xy(nGhost, YY) = y0 + tijGhost[1] * dy_di + tijGhost[2] * dy_dj;

      std::cout << "xy " <<  xy(nGhost, XX)  << " "  << xy(nGhost, YY) << std::endl;

      // Same lonlat as concrete points
      lonlat(nGhost, LON) = lonlat( nOwned, LON);
      lonlat(nGhost, LAT) = lonlat( nOwned, LAT);

      std::cout << "lonlat " <<  lonlat(nGhost, LON)  << " "  <<    lonlat(nGhost, LAT) << std::endl;

      // Is ghost node
      mesh::Nodes::Topology::set(flags(nGhost), mesh::Nodes::Topology::GHOST);
      ghost(nGhost) = 1;

      std::cout << "is ghost node" << std::endl;

      // Partitioning logic to be added in future PR
      glb_idx(nGhost) = nGhost + 1;
      remote_idx(nGhost) = nGhost;

      std::cout << "global remote idx " << std::endl;

      // not sure (below - for multiple PEs)
      part(nGhost) = distribution.partition(nOwned);

      std::cout << "partitioning logic" << std::endl;

      // Append metadata
      // Global indicies of ghost node and owned node
      ghostGblIdx.push_back(nGhost + 1);
      ownedGblIdx.push_back(nOwned + 1);

      std::cout << "append metadata" << std::endl;

      ++nGhost;

      std::cout << "end of loop" << std::endl;

    };

    // Loop over ghost (t, i, j)
    for (auto& pPair : ghostToOwnedTij) addGhostNode(pPair.first, pPair.second);


    // Assert that the correct number of nodes have been set when duplicates are added
    ATLAS_ASSERT( nnodes_all == nGhost, "Insufficient nodes" );

    for ( idx_t it = 0; it < nTiles; it++ ) {
        for ( idx_t ix = 0; ix < N + 1; ix++ ) {
            for ( idx_t iy = 0; iy < N + 1; iy++ ) {
                ATLAS_ASSERT( NodeArray( it, ix, iy ) != -9999, "Node Array Not Set Properly" );
            }
        }
    }

    // Cells in mesh
    mesh.cells().add( new mesh::temporary::Quadrilateral(), nTiles * N * N );
    //int quad_begin  = mesh.cells().elements( 0 ).begin();
    auto cells_part = array::make_view<int, 1>( mesh.cells().partition() );
    auto cells_gidx = array::make_view<gidx_t, 1>( mesh.cells().global_index() );
    auto cells_ridx = array::make_indexview<idx_t, 1>( mesh.cells().remote_index() );
    atlas::mesh::HybridElements::Connectivity& node_connectivity = mesh.cells().node_connectivity();

    int icell = 0;
    idx_t quad_nodes[4];

    for ( int it = 0; it < nTiles; it++ ) {
        for ( int ix = 0; ix < N; ix++ ) {
            for ( int iy = 0; iy < N; iy++ ) {
                quad_nodes[0] = NodeArray( it, ix, iy );
                quad_nodes[1] = NodeArray( it, ix + 1, iy );
                quad_nodes[2] = NodeArray( it, ix + 1, iy + 1 );
                quad_nodes[3] = NodeArray( it, ix, iy + 1 );

                node_connectivity.set( icell, quad_nodes );
                cells_part( icell ) = part( quad_nodes[0] );
                cells_gidx( icell ) = icell + 1;
                cells_ridx( icell ) = icell;

                ++icell;
            }
        }
    }

    // Assertion that correct number of cells are set
    ATLAS_ASSERT( ncells == icell, "Insufficient cells have been set" );

    // Parallel
    generateGlobalElementNumbering( mesh );
    nodes.metadata().set( "parallel", true );

    // Global indices of ghost nodes.
    nodes.metadata().set( "ghost-global-idx", ghostGblIdx);

    // Global indices of owned nodes for each ghost node (same order as above)
    nodes.metadata().set( "owned-global-idx", ownedGblIdx);
}

// -------------------------------------------------------------------------------------------------

void CubedSphereMeshGenerator::hash( eckit::Hash& h ) const {
    h.add( "CubedSphereMeshGenerator" );
    options.hash( h );
}

// -------------------------------------------------------------------------------------------------

namespace {
static MeshGeneratorBuilder<CubedSphereMeshGenerator> CubedSphereMeshGenerator(
    CubedSphereMeshGenerator::static_type() );
}

// -------------------------------------------------------------------------------------------------

}  // namespace meshgenerator
}  // namespace atlas
