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

#include "atlas/array.h"
#include "atlas/field/Field.h"
#include "atlas/grid/CubedSphereGrid.h"
#include "atlas/grid/Distribution.h"
#include "atlas/grid/Partitioner.h"
#include "atlas/library/config.h"
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
      throw_Exception( "CubedSphereMeshGenerator can only work with a cubed-sphere grid", Here() );
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

  int cubeNx = csgrid.GetCubeNx();

  // Number of nodes
  int nnodes_unique = 6*cubeNx*cubeNx+2;        // Number of unique grid nodes
  int nnodes_wghost = 6*(cubeNx+1)*(cubeNx+1);  // With ghost, i.e. including repeat nodes along edges

  // Number of cells
  int ncells = 6*cubeNx*cubeNx;                 // Number of grid cells

  // Construct mesh nodes
  // --------------------
  mesh.nodes().resize( nnodes_wghost );
  mesh::Nodes& nodes = mesh.nodes();
  auto xy         = array::make_view<double, 2>( nodes.xy() );
  auto lonlat     = array::make_view<double, 2>( nodes.lonlat() );
  auto glb_idx    = array::make_view<gidx_t, 1>( nodes.global_index() );
  auto remote_idx = array::make_indexview<idx_t, 1>( nodes.remote_index() );
  auto part       = array::make_view<int, 1>( nodes.partition() );
  auto ghost      = array::make_view<int, 1>( nodes.ghost() );
  auto flags      = array::make_view<int, 1>( nodes.flags() );

  // Start end points for this processor
  // -----------------------------------
  int isc = 0;
  int iec = cubeNx+1;
  int jsc = 0;
  int jec = cubeNx+1;

  std::vector<int> tile_array_me;

  tile_array_me.push_back(0);
  tile_array_me.push_back(1);
  tile_array_me.push_back(2);
  tile_array_me.push_back(3);
  tile_array_me.push_back(4);
  tile_array_me.push_back(5);


  // Array to hold whether a ghost point
  // -----------------------------------
  array::ArrayT<int> isGhostArray( cubeNx+1, cubeNx+1, 6 );
  array::ArrayView<int, 3> isGhost = array::make_view<int, 3>( isGhostArray );

 // assign?

  // Default is not ghost
  for ( int ix = 0; ix < cubeNx+1; ix++ ) {
    for ( int iy = 0; iy < cubeNx+1; iy++ ) {
      for ( int it = 0; it < 6; it++ ) {
        isGhost(ix, iy, it) = 0;
      }
    }
  }

  // Far edges are ghost point
  for ( int it = 0; it < 6; it++ ) {
    for ( int ix = 0; ix < cubeNx+1; ix++ ) {
      isGhost(cubeNx+1, ix, it) = 1;
      isGhost(ix, cubeNx+1, it) = 1;
    }
  }

  // Two special points that are not ghost
  isGhost(0,cubeNx+1,0) = 0; // First tile
  isGhost(1,0,cubeNx+1) = 0; // Second tile

  // Loop over entire grid, including the ghost points
  // -------------------------------------------------

  int inode = 0;

  array::ArrayT<int> NodeArrayT( 6, cubeNx+1, cubeNx+1 );
  array::ArrayView<int, 3> NodeArray = array::make_view<int, 3>( NodeArrayT );

  for ( int it = 0; it < tile_array_me.size(); it++ ) {
    for ( int ix = isc; ix < iec; ix++ ) {
      for ( int iy = jsc; iy < jec; iy++ ) {

        // Get xy from global xy grid array
        double xy_[3];
        csgrid.xy( ix, iy, it, xy_ );

        xy( inode, LON ) = xy_[LON];
        xy( inode, LAT ) = xy_[LAT];

        double lonlat_[2];

        std::cout << "ix, iy, it: " << ix << ", " << iy << ", " << it << std::endl;

        csgrid.lonlat( ix, iy, it, lonlat_ );

        lonlat( inode, LON ) = lonlat_[LON];
        lonlat( inode, LAT ) = lonlat_[LAT];

        // Ghost nodes
        ghost(inode) = isGhost(ix, iy, it);

        glb_idx(inode) = inode;
        remote_idx(inode) = inode;
        part(inode) = 0;

        NodeArray(it, ix, iy) = inode;

        ++inode;

      }
    }
  }

  std::cout << "Begin connectivity" << std::endl;

  // Cells in mesh
  mesh.cells().add( new mesh::temporary::Quadrilateral(), ncells );
  //int quad_begin  = mesh.cells().elements( 0 ).begin();
  auto cells_part = array::make_view<int, 1>( mesh.cells().partition() );
  auto& node_connectivity = mesh.cells().node_connectivity();


  int icell = 0;
  idx_t quad_nodes[4];

  std::cout << "Begin connectivity loop" << std::endl;

  for ( int it = 0; it < tile_array_me.size(); it++ ) {
    for ( int ix = isc; ix < iec-1; ix++ ) {
      for ( int iy = jsc; iy < jec-1; iy++ ) {

        quad_nodes[0] = NodeArray(it, ix  , iy  );
        quad_nodes[1] = NodeArray(it, ix  , iy+1);
        quad_nodes[2] = NodeArray(it, ix+1, iy+1);
        quad_nodes[3] = NodeArray(it, ix+1, iy  );

        node_connectivity.set( icell, quad_nodes );
        cells_part( icell ) = 0;

        ++icell;

      }
    }
  }

  std::cout << "Begin Parallel" << std::endl;

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
  static MeshGeneratorBuilder<CubedSphereMeshGenerator> CubedSphereMeshGenerator( "cubedsphere" );
}

// -------------------------------------------------------------------------------------------------

}  // namespace meshgenerator
}  // namespace atlas
