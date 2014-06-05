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
#include "atlas/meshgen/EqualAreaPartitioner.hpp"
#include "atlas/io/Gmsh.hpp"
#include "atlas/mesh/Mesh.hpp"
#include "atlas/mesh/FunctionSpace.hpp"
#include "atlas/mesh/Field.hpp"
#include "atlas/mesh/Metadata.hpp"
#include "atlas/mesh/Array.hpp"
#include "atlas/mesh/ArrayView.hpp"
#include "atlas/mesh/IndexView.hpp"
#include "atlas/actions/BuildParallelFields.hpp"
#include "atlas/mesh/Parameters.hpp"

using namespace atlas;
using namespace atlas::meshgen;

#define DISABLE if(0)
#define ENABLE if(1)

namespace atlas {
namespace test {

struct IsGhost
{
  IsGhost( FunctionSpace& nodes )
  {
    part    = ArrayView<int,1> (nodes.field("partition") );
    loc_idx = IndexView<int,1> (nodes.field("remote_loc_idx") );
    mypart  = MPL::rank();
  }
  bool operator()(int idx)
  {
    if( part   [idx] != mypart ) return true;
    if( loc_idx[idx] != idx    ) return true;
    return false;
  }
  int mypart;
  ArrayView<int,1> part;
  IndexView<int,1> loc_idx;
};

} // end namespace test
} // end namespace atlas

BOOST_AUTO_TEST_CASE( init ) { MPL::init(); }

BOOST_AUTO_TEST_CASE( test1 )
{
  Mesh* m = new Mesh();
  m->add_function_space( new FunctionSpace( "nodes", "shapefunc", Extents(10,Field::UNDEF_VARS) ) );
  m->function_space("nodes").create_field<double>("coordinates",2);
  actions::build_parallel_fields(*m);

  BOOST_CHECK( m->function_space("nodes").has_field("glb_idx") );
  BOOST_CHECK( m->function_space("nodes").has_field("partition") );
  BOOST_CHECK( m->function_space("nodes").has_field("remote_loc_idx") );

  test::IsGhost is_ghost( m->function_space("nodes") );
  is_ghost.part[0] = MPL::rank();    is_ghost.loc_idx[0] = 0;
  is_ghost.part[1] = MPL::rank();    is_ghost.loc_idx[1] = 0;
  is_ghost.part[2] = MPL::rank()+1;  is_ghost.loc_idx[2] = 2;

  BOOST_CHECK_EQUAL( is_ghost(0), false );
  BOOST_CHECK_EQUAL( is_ghost(1), true );
  BOOST_CHECK_EQUAL( is_ghost(2), true );

  //std::stringstream filename; filename << "T63_halo_EW.msh";
  //Gmsh::write(*m,filename.str());
  delete m;
}

BOOST_AUTO_TEST_CASE( test2 )
{
  RGGMeshGenerator generate;
  Mesh* m = generate( T63() );

  actions::build_parallel_fields(*m);

  FunctionSpace& nodes = m->function_space("nodes");
  IndexView<int,1> loc_idx ( nodes.field("remote_loc_idx") );
  ArrayView<int,1> part    ( nodes.field("partition")      );
  ArrayView<int,1> glb_idx ( nodes.field("glb_idx")        );

  test::IsGhost is_ghost(nodes);

  int nb_periodic = 0;
  for( int jnode=0; jnode<nodes.extents()[0]; ++jnode )
  {

//    std::cout << part( jnode ) <<"["<<loc_idx( jnode ) << "] master of " << MPL::rank() << "["<<jnode<<"]" << std::endl;

    if( is_ghost(jnode) ) ++nb_periodic;
  }

  BOOST_CHECK_EQUAL( nb_periodic, 64 );

  std::stringstream filename; filename << "periodic.msh";
  Gmsh::write(*m,filename.str());
  delete m;
}

BOOST_AUTO_TEST_CASE( finalize ) { MPL::finalize(); }
