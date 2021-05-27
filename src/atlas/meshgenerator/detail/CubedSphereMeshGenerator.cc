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
#include "atlas/meshgenerator/detail/MeshGeneratorFactory.h"
#include "atlas/meshgenerator/detail/CubedSphereMeshGenerator.h"
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

void CubedSphereMeshGenerator::generate( const Grid& grid, const grid::Distribution& distribution,
                                         Mesh& mesh ) const {

  const auto csgrid = CubedSphereGrid( grid );

  const int N = csgrid.N();
  const int nTiles = csgrid.GetNTiles();

  ATLAS_TRACE("CubedSphereMeshGenerator::generate");
  Log::debug() << "Number of faces per tile edge = " << std::to_string(N) << std::endl;

  // Number of nodes
  int nnodes = nTiles*N*N+2;             // Number of unique grid nodes
  int nnodes_all = nTiles*(N+1)*(N+1);   // Number of grid nodes including edge and corner duplicates
  int ncells = nTiles*N*N;               // Number of unique grid cells

  // Construct mesh nodes
  // --------------------
  mesh.nodes().resize( nnodes );
  mesh::Nodes& nodes = mesh.nodes();
  auto xy         = array::make_view<double, 2>( nodes.xy() );
  auto lonlat     = array::make_view<double, 2>( nodes.lonlat() );
  auto glb_idx    = array::make_view<gidx_t, 1>( nodes.global_index() );
  auto remote_idx = array::make_indexview<idx_t, 1>( nodes.remote_index() );
  auto part       = array::make_view<int, 1>( nodes.partition() );
  auto ghost      = array::make_view<int, 1>( nodes.ghost() );
  //auto flags      = array::make_view<int, 1>( nodes.flags() );

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

  array::ArrayT<int> NodeArrayT( nTiles, N+1, N+1 ); // All grid points including duplicates
  auto NodeArray = array::make_view<int, 3>( NodeArrayT );


  {
      ATLAS_TRACE("iterate entire grid");
      idx_t n{0};
      for( auto& p : csgrid.tij() ) {
          NodeArray(p.t(), p.i(), p.j()) = n++;
      }
  }


  for ( it = 0; it < nTiles; it++ ) {               // 0, 1, 2, 3, 4, 5
    for ( ix = 0; ix < N; ix++ ) {        // 0, 1, ..., N-1
      for ( iy = 0; iy < N; iy++ ) {      // 0, 1, ..., N-1

        // Get xy from global xy grid array
        csgrid.xy( ix, iy, it, xy_ );

        idx_t n = NodeArray(it, ix, iy);

        xy( n, XX ) = xy_[XX];
        xy( n, YY ) = xy_[YY];

        csgrid.lonlat( ix, iy, it, lonlat_ );

        lonlat( n, LON ) = lonlat_[LON];
        lonlat( n, LAT ) = lonlat_[LAT];

        // Ghost nodes
        ghost(n) = 0; // No ghost nodes

        glb_idx(n) = n+1;
        remote_idx(n) = n;
        part(n) = distribution.partition(n);
        inode++;
      }
    }
  }

  // Extra point 1
  // -------------
  {
      it = 0;
      ix = 0;
      iy = N;

      idx_t n = NodeArray(it, ix, iy);
      csgrid.xy( ix, iy, it, xy_ );
      xy( n, XX ) = xy_[XX];
      xy( n, YY ) = xy_[YY];
      csgrid.lonlat( ix, iy, it, lonlat_ );
      lonlat( n, LON ) = lonlat_[LON];
      lonlat( n, LAT ) = lonlat_[LAT];
      ghost(n) = 0;
      glb_idx(n) = n+1;
      remote_idx(n) = n;
      part(n) = distribution.partition(n);
      inode++;
  }

  // Extra point 2
  // -------------
  {
      it = 1;
      ix = N;
      iy = 0;
      idx_t n = NodeArray(it, ix, iy);

      csgrid.xy( ix, iy, it, xy_ );
      xy( n, XX ) = xy_[XX];
      xy( n, YY ) = xy_[YY];
      csgrid.lonlat( ix, iy, it, lonlat_ );
      lonlat( n, LON ) = lonlat_[LON];
      lonlat( n, LAT ) = lonlat_[LAT];
      ghost(n) = 0;
      glb_idx(n) = n+1;
      remote_idx(n) = n;
      part(n) = distribution.partition(n);
      inode++;
  }

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

  // Fill duplicate points in the corners
  // ------------------------------------
  NodeArray(0, N, N) = NodeArray(2, 0, 0); ++inode;
  NodeArray(1, N, N) = NodeArray(3, 0, 0); ++inode;
  NodeArray(2, N, N) = NodeArray(4, 0, 0); ++inode;
  NodeArray(3, N, N) = NodeArray(5, 0, 0); ++inode;
  NodeArray(4, N, N) = NodeArray(0, 0, 0); ++inode;
  NodeArray(5, N, N) = NodeArray(1, 0, 0); ++inode;

  // Special points have two duplicates each
  NodeArray(2, 0, N) = NodeArray(0, 0, N); ++inode;
  NodeArray(4, 0, N) = NodeArray(0, 0, N); ++inode;
  NodeArray(3, N, 0) = NodeArray(1, N, 0); ++inode;
  NodeArray(5, N, 0) = NodeArray(1, N, 0); ++inode;

  // Top & right duplicates
  // ----------------------

  // Tile 1
  for ( ix = 1; ix < N; ix++ ) {
    NodeArray(0, ix, N) = NodeArray(2, 0, N-ix); ++inode;
  }
  for ( iy = 0; iy < N; iy++ ) {
    NodeArray(0, N, iy) = NodeArray(1, 0, iy); ++inode;
  }

  // Tile 2
  for ( ix = 0; ix < N; ix++ ) {
    NodeArray(1, ix, N) = NodeArray(2, ix, 0); ++inode;
  }
  for ( iy = 1; iy < N; iy++ ) {
    NodeArray(1, N, iy) = NodeArray(3, N-iy, 0); ++inode;
  }

  // Tile 3
  for ( ix = 1; ix < N; ix++ ) {
    NodeArray(2, ix, N) = NodeArray(4, 0, N-ix); ++inode;
  }
  for ( iy = 0; iy < N; iy++ ) {
    NodeArray(2, N, iy) = NodeArray(3, 0, iy); ++inode;
  }

  // Tile 4
  for ( ix = 0; ix < N; ix++ ) {
    NodeArray(3, ix, N) = NodeArray(4, ix, 0); ++inode;
  }
  for ( iy = 1; iy < N; iy++ ) {
    NodeArray(3, N, iy) = NodeArray(5, N-iy, 0); ++inode;
  }

  // Tile 5
  for ( ix = 1; ix < N; ix++ ) {
    NodeArray(4, ix, N) = NodeArray(0, 0, N-ix); ++inode;
  }
  for ( iy = 0; iy < N; iy++ ) {
    NodeArray(4, N, iy) = NodeArray(5, 0, iy); ++inode;
  }

  // Tile 6
  for ( ix = 0; ix < N; ix++ ) {
    NodeArray(5, ix, N) = NodeArray(0, ix, 0); ++inode;
  }
  for ( iy = 1; iy < N; iy++ ) {
    NodeArray(5, N, iy) = NodeArray(1, N-iy, 0); ++inode;
  }

  // Assert that the correct number of nodes have been set when duplicates are added
  ATLAS_ASSERT( nnodes_all == inode, "Insufficient nodes" );

  for ( it = 0; it < nTiles; it++ ) {
    for ( ix = 0; ix < N+1; ix++ ) {
      for ( iy = 0; iy < N+1; iy++ ) {
        ATLAS_ASSERT( NodeArray(it, ix, iy) != -9999, "Node Array Not Set Properly" );
      }
    }
  }

  // Cells in mesh
  mesh.cells().add( new mesh::temporary::Quadrilateral(), nTiles*N*N );
  //int quad_begin  = mesh.cells().elements( 0 ).begin();
  auto cells_part    = array::make_view<int, 1>( mesh.cells().partition() );
  auto cells_gidx = array::make_view<gidx_t, 1>( mesh.cells().global_index() );
  auto cells_ridx = array::make_indexview<idx_t, 1>( mesh.cells().remote_index() );
  atlas::mesh::HybridElements::Connectivity& node_connectivity = mesh.cells().node_connectivity();

  int icell = 0;
  idx_t quad_nodes[4];

  for ( int it = 0; it < nTiles; it++ ) {
    for ( int ix = 0; ix < N; ix++ ) {
      for ( int iy = 0; iy < N; iy++ ) {

        quad_nodes[0] = NodeArray(it, ix  , iy  );
        quad_nodes[1] = NodeArray(it, ix+1, iy  );
        quad_nodes[2] = NodeArray(it, ix+1, iy+1);
        quad_nodes[3] = NodeArray(it, ix  , iy+1);

        node_connectivity.set( icell, quad_nodes );
        cells_part( icell ) = part( quad_nodes[0] );
        cells_gidx( icell ) = icell + 1 ;
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
  static MeshGeneratorBuilder<CubedSphereMeshGenerator> CubedSphereMeshGenerator( CubedSphereMeshGenerator::static_type() );
}

// -------------------------------------------------------------------------------------------------

}  // namespace meshgenerator
}  // namespace atlas
