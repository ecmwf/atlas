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
#define BOOST_UNIT_TEST_FRAMEWORK_HEADER_ONLY
#include "ecbuild/boost_test_framework.h"

#include "atlas/mpl/MPL.hpp"
#include "atlas/atlas_config.h"
#include "atlas/meshgen/RGG.hpp"
#include "atlas/meshgen/RGGMeshGenerator.hpp"
#include "atlas/meshgen/EqualAreaPartitioner.hpp"
#include "atlas/io/Gmsh.hpp"
#include "atlas/mesh/Mesh.hpp"
#include "atlas/mesh/FunctionSpace.hpp"
#include "atlas/mesh/Field.hpp"
#include "atlas/mesh/Metadata.hpp"
#include "atlas/util/Array.hpp"
#include "atlas/util/ArrayView.hpp"
#include "atlas/util/IndexView.hpp"
#include "atlas/actions/BuildParallelFields.hpp"
#include "atlas/actions/BuildPeriodicBoundaries.hpp"
#include "atlas/mesh/Parameters.hpp"
#include "atlas/mesh/Util.hpp"

using namespace atlas;
using namespace atlas::meshgen;

#define DISABLE if(0)
#define ENABLE if(1)

BOOST_AUTO_TEST_CASE( init ) { MPL::init(); }

BOOST_AUTO_TEST_CASE( test1 )
{
  Mesh* m = new Mesh();
  FunctionSpace* nodes =  new FunctionSpace( "nodes", "shapefunc", Extents(10,Field::UNDEF_VARS) );
  m->add_function_space( nodes );
  ArrayView<double,2> latlon( m->function_space("nodes").create_field<double>("coordinates",2) );
  ArrayView<int,1> glb_idx( m->function_space("nodes").create_field<int>("glb_idx",1) );
  ArrayView<int,1> part ( nodes->create_field<int>("partition",1) );

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

  latlon(0,XX) = 0.;         latlon(0,YY) = 1.;
  latlon(1,XX) = 0.;         latlon(1,YY) =-1.;
  latlon(2,XX) = M_PI/4.;    latlon(2,YY) = 1.;
  latlon(3,XX) = M_PI/4.;    latlon(3,YY) =-1.;
  latlon(4,XX) = M_PI/2.;    latlon(4,YY) = 1.;
  latlon(5,XX) = M_PI/2.;    latlon(5,YY) =-1.;
  latlon(6,XX) = 3.*M_PI/4.; latlon(6,YY) = 1.;
  latlon(7,XX) = 3.*M_PI/4.; latlon(7,YY) =-1.;
  latlon(8,XX) = 2.*M_PI;    latlon(8,YY) = 1.;
  latlon(9,XX) = 2.*M_PI;    latlon(9,YY) =-1.;

  actions::build_parallel_fields(*m);

  BOOST_REQUIRE( m->function_space("nodes").has_field("remote_idx") );

  IndexView<int,1> loc( nodes->field("remote_idx") );
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

  switch ( MPL::rank() )
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
  if( MPL::rank() == 1 )
  {
    BOOST_CHECK_EQUAL( is_ghost(8), true );
    BOOST_CHECK_EQUAL( is_ghost(9), true );
  }

  delete m;
}

BOOST_AUTO_TEST_CASE( test2 )
{
  RGGMeshGenerator generate;
  generate.options.set("nb_parts",MPL::size());
  generate.options.set("part",MPL::rank());
  Mesh* m = generate( T63() );
  actions::build_parallel_fields(*m);

  FunctionSpace& nodes = m->function_space("nodes");
  IndexView<int,1> loc_idx ( nodes.field("remote_idx") );
  ArrayView<int,1> part    ( nodes.field("partition")      );
  ArrayView<int,1> glb_idx ( nodes.field("glb_idx")        );

  IsGhost is_ghost(nodes);

  int nb_ghost = 0;
  for( int jnode=0; jnode<nodes.extents()[0]; ++jnode )
  {
    if( is_ghost(jnode) ) ++nb_ghost;
  }

  if( MPL::rank() == 0 ) BOOST_CHECK_EQUAL( nb_ghost, 129 );
  if( MPL::rank() == 1 ) BOOST_CHECK_EQUAL( nb_ghost, 0   );

  actions::build_periodic_boundaries(*m);

  int nb_periodic = -nb_ghost;
  for( int jnode=0; jnode<nodes.extents()[0]; ++jnode )
  {
    if( is_ghost(jnode) ) ++nb_periodic;
  }

  if( MPL::rank() == 0 ) BOOST_CHECK_EQUAL( nb_periodic, 32 );
  if( MPL::rank() == 1 ) BOOST_CHECK_EQUAL( nb_periodic, 32 );

  std::stringstream filename; filename << "periodic_p"<<MPL::rank()<<".msh";
  Gmsh().write(*m,filename.str());
  delete m;
}


BOOST_AUTO_TEST_CASE( finalize ) { MPL::finalize(); }
