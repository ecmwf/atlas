/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <sstream>
#include <algorithm>

#define BOOST_TEST_MODULE TestBuildParallelFields
#include "ecbuild/boost_test_framework.h"

#include "eckit/config/ResourceMgr.h"
#include "atlas/actions/BuildParallelFields.h"
#include "atlas/actions/BuildPeriodicBoundaries.h"
#include "atlas/atlas_config.h"
#include "atlas/io/Gmsh.h"
#include "atlas/Field.h"
#include "atlas/FunctionSpace.h"
#include "atlas/Mesh.h"
#include "atlas/Metadata.h"
#include "atlas/Parameters.h"
#include "atlas/Util.h"
#include "atlas/meshgen/EqualRegionsPartitioner.h"
#include "atlas/grids/grids.h"
#include "atlas/meshgen/ReducedGridMeshGenerator.h"
#include "atlas/mpi/mpi.h"
#include "atlas/util/Array.h"
#include "atlas/util/ArrayView.h"
#include "atlas/util/IndexView.h"

using namespace atlas::io;
using namespace atlas::meshgen;

namespace atlas {
namespace test {

#define DISABLE if(0)
#define ENABLE if(1)

BOOST_AUTO_TEST_CASE( init )
{
  eckit::mpi::init();
  eckit::ResourceMgr::instance().set("atlas.meshgen.angle","27.5");
  eckit::ResourceMgr::instance().set("atlas.meshgen.triangulate","false");
}

BOOST_AUTO_TEST_CASE( test1 )
{
	Mesh::Ptr m = Mesh::create();

	FunctionSpace& nodes = m->create_function_space( "nodes", "shapefunc", make_shape(10,Field::UNDEF_VARS) );

  ArrayView<double,2> lonlat ( m->function_space("nodes").create_field<double>("lonlat",2) );
  ArrayView<gidx_t,1> glb_idx( m->function_space("nodes").create_field<gidx_t>("glb_idx",1) );
  ArrayView<int,1> part ( nodes.create_field<int>("partition",1) );
  ArrayView<int,1> flags ( nodes.create_field<int>("flags",1) );
  flags = Topology::NONE;

  // This is typically available
  glb_idx(0) = 1;    part(0) = 0;
  glb_idx(1) = 2;    part(1) = 0;
  glb_idx(2) = 3;    part(2) = 0;
  glb_idx(3) = 4;    part(3) = 0;
  glb_idx(4) = 5;    part(4) = 0;
  glb_idx(5) = 6;    part(5) = 1;
  glb_idx(6) = 7;    part(6) = 1;
  glb_idx(7) = 8;    part(7) = 1;
  glb_idx(8) = 9;    part(8) = 1;
  glb_idx(9) = 10;   part(9) = 1;

  lonlat(0,LON) = 0.;    lonlat(0,LAT) = 80.;    Topology::set( flags(0), Topology::BC|Topology::WEST );
  lonlat(1,LON) = 0.;    lonlat(1,LAT) =-80.;    Topology::set( flags(1), Topology::BC|Topology::WEST );
  lonlat(2,LON) = 90.;   lonlat(2,LAT) = 80.;
  lonlat(3,LON) = 90.;   lonlat(3,LAT) =-80.;
  lonlat(4,LON) = 180.;  lonlat(4,LAT) = 80.;
  lonlat(5,LON) = 180.;  lonlat(5,LAT) =-80.;
  lonlat(6,LON) = 270.;  lonlat(6,LAT) = 80.;
  lonlat(7,LON) = 270.;  lonlat(7,LAT) =-80.;
  lonlat(8,LON) = 360.;  lonlat(8,LAT) = 80.;    Topology::set( flags(8), Topology::BC|Topology::EAST );
  lonlat(9,LON) = 360.;  lonlat(9,LAT) =-80.;    Topology::set( flags(9), Topology::BC|Topology::EAST );

  actions::build_parallel_fields(*m);

  BOOST_REQUIRE( m->function_space("nodes").has_field("remote_idx") );

  IndexView<int,1> loc( nodes.field("remote_idx") );
  BOOST_CHECK_EQUAL( loc(0) , 0 );
  BOOST_CHECK_EQUAL( loc(1) , 1 );
  BOOST_CHECK_EQUAL( loc(2) , 2 );
  BOOST_CHECK_EQUAL( loc(3) , 3 );
  BOOST_CHECK_EQUAL( loc(4) , 4 );
  BOOST_CHECK_EQUAL( loc(5) , 5 );
  BOOST_CHECK_EQUAL( loc(6) , 6 );
  BOOST_CHECK_EQUAL( loc(7) , 7 );
  BOOST_CHECK_EQUAL( loc(8) , 8 );
  BOOST_CHECK_EQUAL( loc(9) , 9 );

  IsGhost is_ghost( m->function_space("nodes") );

  switch ( eckit::mpi::rank() )
  {
  case 0:
    BOOST_CHECK_EQUAL( is_ghost(0), false );
    BOOST_CHECK_EQUAL( is_ghost(1), false );
    BOOST_CHECK_EQUAL( is_ghost(2), false );
    BOOST_CHECK_EQUAL( is_ghost(3), false );
    BOOST_CHECK_EQUAL( is_ghost(4), false );
    BOOST_CHECK_EQUAL( is_ghost(5), true );
    BOOST_CHECK_EQUAL( is_ghost(6), true );
    BOOST_CHECK_EQUAL( is_ghost(7), true );
    BOOST_CHECK_EQUAL( is_ghost(8), true );
    BOOST_CHECK_EQUAL( is_ghost(9), true );
    break;
  case 1:
    BOOST_CHECK_EQUAL( is_ghost(0), true );
    BOOST_CHECK_EQUAL( is_ghost(1), true );
    BOOST_CHECK_EQUAL( is_ghost(2), true );
    BOOST_CHECK_EQUAL( is_ghost(3), true );
    BOOST_CHECK_EQUAL( is_ghost(4), true );
    BOOST_CHECK_EQUAL( is_ghost(5), false );
    BOOST_CHECK_EQUAL( is_ghost(6), false );
    BOOST_CHECK_EQUAL( is_ghost(7), false );
    BOOST_CHECK_EQUAL( is_ghost(8), false );
    BOOST_CHECK_EQUAL( is_ghost(9), false );
    break;
  }

  actions::build_periodic_boundaries(*m);

  // points 8 and 9 are periodic slave points of points 0 and 1
  BOOST_CHECK_EQUAL( part(8), 0 );
  BOOST_CHECK_EQUAL( part(9), 0 );
  BOOST_CHECK_EQUAL( loc(8), 0 );
  BOOST_CHECK_EQUAL( loc(9), 1 );
  if( eckit::mpi::rank() == 1 )
  {
    BOOST_CHECK_EQUAL( is_ghost(8), true );
    BOOST_CHECK_EQUAL( is_ghost(9), true );
  }
}

BOOST_AUTO_TEST_CASE( test2 )
{
  ReducedGridMeshGenerator generate;
  generate.options.set<int>("nb_parts",eckit::mpi::size());
  generate.options.set<int>("part",eckit::mpi::rank());
  Mesh* m = generate( grids::rgg::N32() );
  actions::build_parallel_fields(*m);

  FunctionSpace& nodes = m->function_space("nodes");
  IndexView<int,1> loc_idx ( nodes.field("remote_idx") );
  ArrayView<int,1> part    ( nodes.field("partition")      );
  ArrayView<gidx_t,1> glb_idx ( nodes.field("glb_idx")        );

  IsGhost is_ghost(nodes);

  int nb_ghost = 0;
  for( int jnode=0; jnode<nodes.shape(0); ++jnode )
  {
    if( is_ghost(jnode) ) ++nb_ghost;
  }

  if( eckit::mpi::rank() == 0 ) BOOST_CHECK_EQUAL( nb_ghost, 129 );
  if( eckit::mpi::rank() == 1 ) BOOST_CHECK_EQUAL( nb_ghost, 0   );

  actions::build_periodic_boundaries(*m);

  int nb_periodic = -nb_ghost;
  for( int jnode=0; jnode<nodes.shape(0); ++jnode )
  {
    if( is_ghost(jnode) ) ++nb_periodic;
  }

  if( eckit::mpi::rank() == 0 ) BOOST_CHECK_EQUAL( nb_periodic, 32 );
  if( eckit::mpi::rank() == 1 ) BOOST_CHECK_EQUAL( nb_periodic, 32 );

  std::stringstream filename; filename << "periodic.msh";
  Gmsh().write(*m,filename.str());
  delete m;
}


BOOST_AUTO_TEST_CASE( finalize ) { eckit::mpi::finalize(); }

} // namespace test
} // namespace atlas
