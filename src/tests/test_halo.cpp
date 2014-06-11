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

#define BOOST_TEST_MODULE TestHalo
#define BOOST_UNIT_TEST_FRAMEWORK_HEADER_ONLY
#include "ecbuild/boost_test_framework.h"

#include "atlas/mpl/MPL.hpp"
#include "atlas/atlas_config.h"
#include "tests/TestMeshes.hpp"
#include "atlas/io/Gmsh.hpp"
#include "atlas/mesh/Mesh.hpp"
#include "atlas/mesh/FunctionSpace.hpp"
#include "atlas/mesh/IndexView.hpp"
#include "atlas/mesh/ArrayView.hpp"
#include "atlas/actions/BuildHalo.hpp"
#include "atlas/actions/BuildParallelFields.hpp"
#include "atlas/mesh/Parameters.hpp"

using namespace atlas;
using namespace atlas::meshgen;

struct MPIFixture {
    MPIFixture()  { MPL::init(); }
    ~MPIFixture() { MPL::finalize(); }
};

BOOST_GLOBAL_FIXTURE( MPIFixture )

BOOST_AUTO_TEST_CASE( test_small )
{
  int nlat = 5;
  int lon[5] = {10, 12, 14, 16, 16};
  Mesh::Ptr m = test::generate_mesh(nlat, lon);

  actions::build_parallel_fields(*m);
  actions::make_periodic(*m);
  actions::build_halo(*m,2);

  IndexView<int,1> ridx ( m->function_space("nodes").field("remote_idx") );
  ArrayView<int,1> gidx ( m->function_space("nodes").field("glb_idx") );

  if( MPL::size() == 5 )
  {
    switch( MPL::rank() ) // with 5 tasks
    {
    case 0:
      BOOST_CHECK_EQUAL( ridx(9),  9  );
      BOOST_CHECK_EQUAL( gidx(9),  10 );
      BOOST_CHECK_EQUAL( ridx(29), 9 );
      BOOST_CHECK_EQUAL( gidx(29), 875430066 ); // hashed unique idx
      break;
    }
  }
  else
  {
    if( MPL::rank() == 0 )
      std::cout << "skipping tests with 5 mpi tasks!" << std::endl;
  }

  std::stringstream filename; filename << "small_halo_p" << MPL::rank() << ".msh";
  Gmsh::write(*m,filename.str());
}

#if 1
BOOST_AUTO_TEST_CASE( test_t63 )
{
  Mesh::Ptr m = test::generate_mesh( T63() );

  actions::build_parallel_fields(*m);
  actions::make_periodic(*m);
  actions::build_halo(*m,5);

  std::stringstream filename; filename << "T63_halo_p" << MPL::rank() << ".msh";
  Gmsh::write(*m,filename.str());
}
#endif

