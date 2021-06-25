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

    int inode = 0;
    double xy_[3];
    double lonlat_[2];

    int it;
    int ix;
    int iy;

    // Loop over entire grid
    // ---------------------

    // Node array will include duplicates of the grid points to make it easier to fill up the
    // neighbours. However the mesh will not contain these points so no Ghost nodes are required.
    // We could include the duplicate points in the array if we make them Ghost nodes but it is not
    // clear if this provides any benefit.

    array::ArrayT<int> NodeArrayT( nTiles, N + 1, N + 1 );  // All grid points including duplicates
    auto NodeArray = array::make_view<int, 3>( NodeArrayT );


    {
        ATLAS_TRACE( "iterate entire grid" );
        idx_t n{0};
        for ( auto& p : csgrid.tij() ) {
            NodeArray( p.t(), p.i(), p.j() ) = n++;
        }
    }

    auto addOwnedNode = [&](idx_t it, idx_t ix, idx_t iy) {

        // Get xy from global xy grid array

        double xy_[3];
        csgrid.xy( ix, iy, it, xy_ );

        idx_t iGridIt = NodeArray(it, ix, iy);

        xy( iGridIt, XX ) = xy_[XX];
        xy( iGridIt, YY ) = xy_[YY];

        // Get lonlat from global lonlat array
        double lonlat_[2];
        csgrid.lonlat( ix, iy, it, lonlat_ );

        lonlat( iGridIt, LON ) = lonlat_[LON];
        lonlat( iGridIt, LAT ) = lonlat_[LAT];

        // Is not ghost node
        mesh::Nodes::Topology::reset(flags(iGridIt));
        ghost(iGridIt) = 0;

        glb_idx(iGridIt) = iGridIt + 1;
        remote_idx(iGridIt) = iGridIt;
        part(iGridIt) = distribution.partition(iGridIt);
        inode++;

        return;
    };


    for ( it = 0; it < nTiles; it++ ) {     // 0, 1, 2, 3, 4, 5
        for ( ix = 0; ix < N; ix++ ) {      // 0, 1, ..., N-1
            for ( iy = 0; iy < N; iy++ ) {  // 0, 1, ..., N-1

                // Get xy from global xy grid array
                addOwnedNode(it, ix, iy);
            }
        }
    }

    // Extra point 1
    // -------------
    addOwnedNode(0, 0, N);

    // Extra point 2
    // -------------
    addOwnedNode(1, N, 0);

#if 0
  // TESTING: Check that inverse projection returns the correct xyt
  // --------------------------------------------------------------
  idx_t ijtnode[3];
  idx_t ijtproj[3];
  double lonlattest[2];

  for (int n = 0; n < inode-2; n++) {
    ijtnode[0] = -1;
    ijtnode[1] = -1;
    ijtnode[2] = -1;
    ijtproj[0] = -2;
    ijtproj[1] = -2;
    ijtproj[2] = -2;
    lonlattest[LON] = lonlat( n, LON );
    lonlattest[LAT] = lonlat( n, LAT );
    // Search the array for this node
    for ( it = 0; it < nTiles; it++ ) {
      for ( ix = 0; ix < N; ix++ ) {
        for ( iy = 0; iy < N; iy++ ) {
          if (NodeArray(it, ix, iy) == n) {
            ijtnode[0] = ix;
            ijtnode[1] = iy;
            ijtnode[2] = it;
            csgrid.lonlat2xy( lonlattest, ijtproj);
          }
        }
      }
    }
    bool samePoint = ijtproj[0] == ijtnode[0] && ijtproj[1] == ijtnode[1] && ijtproj[2] == ijtnode[2];
    ATLAS_ASSERT(samePoint, "NOT THE SAME POINT");
  }
#endif

    // Assert that the correct number of nodes have been set
    ATLAS_ASSERT( nnodes == inode, "Insufficient nodes" );

    // Start adding Ghost nodes.

    auto addGhostNode = [&](idx_t itGhost, idx_t ixGhost, idx_t iyGhost,
      idx_t itOwned, idx_t ixOwned, idx_t iyOwned) {

      // Get concrete node id
      auto nOwned = NodeArray(itOwned, ixOwned, iyOwned);

      // Add ghost node to NodeArray
      NodeArray(itGhost, ixGhost, iyGhost) = inode;

      // "Create" ghost xy coordinate.

      // Get Jacobian of coords rtw indicesd
      auto x0 = xy(NodeArray(itGhost, 0, 0), XX);
      auto y0 = xy(NodeArray(itGhost, 0, 0), YY);

      auto dx_di = xy(NodeArray(itGhost, 1, 0), XX) - x0;
      auto dx_dj = xy(NodeArray(itGhost, 0, 1), XX) - x0;
      auto dy_di = xy(NodeArray(itGhost, 1, 0), YY) - y0;
      auto dy_dj = xy(NodeArray(itGhost, 0, 1), YY) - y0;

      // Set xy coordinates
      xy(inode, XX) = x0 + ixGhost * dx_di + iyGhost * dx_dj;
      xy(inode, YY) = y0 + ixGhost * dy_di + iyGhost * dy_dj;

      // Same lonlat as concrete points
      lonlat(inode, LON) = lonlat( nOwned, LON);
      lonlat(inode, LAT) = lonlat( nOwned, LAT);

      // Is ghost node
      mesh::Nodes::Topology::set(flags(inode), mesh::Nodes::Topology::GHOST);
      ghost(inode) = 1;

      // Partitioning logic to be added in future PR
      glb_idx(inode) = inode + 1;
      remote_idx(inode) = inode;
      part(inode) = distribution.partition(inode);

      inode++;

    };

    // Fill duplicate points in the corners
    // ------------------------------------
    addGhostNode(0, N, N, 2, 0, 0);
    addGhostNode(1, N, N, 3, 0, 0);
    addGhostNode(2, N, N, 4, 0, 0);
    addGhostNode(3, N, N, 5, 0, 0);
    addGhostNode(4, N, N, 0, 0, 0);
    addGhostNode(5, N, N, 1, 0, 0);

    // Special points have two duplicates each
    addGhostNode(2, 0, N, 0, 0, N);
    addGhostNode(4, 0, N, 0, 0, N);
    addGhostNode(3, N, 0, 1, N, 0);
    addGhostNode(5, N, 0, 1, N, 0);

    // Top & right duplicates
    // ----------------------

    // Tile 1
    for ( ix = 1; ix < N; ix++ ) {
      addGhostNode(0, ix, N, 2, 0, N - ix);
    }
    for ( iy = 0; iy < N; iy++ ) {
      addGhostNode(0, N, iy, 1, 0, iy);
    }

    // Tile 2
    for ( ix = 0; ix < N; ix++ ) {
      addGhostNode(1, ix, N, 2, ix, 0);
    }
    for ( iy = 1; iy < N; iy++ ) {
      addGhostNode(1, N, iy, 3, N - iy, 0);
    }

    // Tile 3
    for ( ix = 1; ix < N; ix++ ) {
      addGhostNode(2, ix, N, 4, 0, N - ix);
    }
    for ( iy = 0; iy < N; iy++ ) {
      addGhostNode(2, N, iy, 3, 0, iy);
    }

    // Tile 4
    for ( ix = 0; ix < N; ix++ ) {
      addGhostNode(3, ix, N, 4, ix, 0);
    }
    for ( iy = 1; iy < N; iy++ ) {
      addGhostNode(3, N, iy, 5, N - iy, 0);
    }

    // Tile 5
    for ( ix = 1; ix < N; ix++ ) {
      addGhostNode(4, ix, N, 0, 0, N - ix);
    }
    for ( iy = 0; iy < N; iy++ ) {
      addGhostNode(4, N, iy, 5, 0, iy);
    }

    // Tile 6
    for ( ix = 0; ix < N; ix++ ) {
      addGhostNode(5, ix, N, 0, ix, 0);
    }
    for ( iy = 1; iy < N; iy++ ) {
      addGhostNode(5, N, iy, 1, N - iy, 0);
    }

    // Assert that the correct number of nodes have been set when duplicates are added
    ATLAS_ASSERT( nnodes_all == inode, "Insufficient nodes" );

    for ( it = 0; it < nTiles; it++ ) {
        for ( ix = 0; ix < N + 1; ix++ ) {
            for ( iy = 0; iy < N + 1; iy++ ) {
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
