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
#include <iomanip>

#define BOOST_TEST_MODULE TestHalo
#define BOOST_UNIT_TEST_FRAMEWORK_HEADER_ONLY
#include "ecbuild/boost_test_framework.h"
#include "eckit/filesystem/LocalPathName.h"

#include "atlas/mpl/MPL.hpp"
#include "atlas/atlas_config.h"
#include "tests/TestMeshes.hpp"
#include "atlas/io/Gmsh.hpp"
#include "atlas/mesh/Mesh.hpp"
#include "atlas/mesh/FunctionSpace.hpp"
#include "atlas/util/IndexView.hpp"
#include "atlas/util/ArrayView.hpp"
#include "atlas/util/Array.hpp"
#include "atlas/actions/BuildHalo.hpp"
#include "atlas/actions/BuildParallelFields.hpp"
#include "atlas/actions/BuildPeriodicBoundaries.hpp"
#include "atlas/actions/BuildEdges.hpp"
#include "atlas/actions/BuildDualMesh.hpp"
#include "atlas/mesh/Parameters.hpp"
#include "atlas/mesh/Util.hpp"

using namespace atlas;
using namespace atlas::meshgen;

namespace atlas {
namespace test {

double dual_volume(Mesh& mesh)
{
  FunctionSpace& nodes = mesh.function_space("nodes");
  IsGhost is_ghost_node(nodes);
  int nb_nodes = nodes.extents()[0];
  ArrayView<double,1> dual_volumes ( nodes.field("dual_volumes") );
  ArrayView<int,1> glb_idx ( nodes.field("glb_idx") );
  double area=0;
  for( int node=0; node<nb_nodes; ++node )
  {
    if( ! is_ghost_node(node) )
    {
      area += dual_volumes(node);
    }
  }
  MPL_CHECK_RESULT( MPI_Allreduce( MPI_IN_PLACE, &area, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD ) );
  return area;
}

}
}

struct MPIFixture {
    MPIFixture()  { MPL::init(); }
    ~MPIFixture() { MPL::finalize(); }
};

BOOST_GLOBAL_FIXTURE( MPIFixture )

#if 0
BOOST_AUTO_TEST_CASE( test_small )
{
  int nlat = 5;
  int lon[5] = {10, 12, 14, 16, 16};

  Mesh::Ptr m = test::generate_mesh(nlat, lon);

  actions::build_parallel_fields(*m);
  actions::build_periodic_boundaries(*m);
  actions::build_halo(*m,2);


  if( MPL::size() == 5 )
  {
    IndexView<int,1> ridx ( m->function_space("nodes").field("remote_idx") );
    ArrayView<int,1> gidx ( m->function_space("nodes").field("glb_idx") );

    switch( MPL::rank() ) // with 5 tasks
    {
    case 0:
      BOOST_CHECK_EQUAL( ridx(9),  9  );
      BOOST_CHECK_EQUAL( gidx(9),  10 );
      BOOST_CHECK_EQUAL( ridx(30), 9 );
      BOOST_CHECK_EQUAL( gidx(30), 875430066 ); // hashed unique idx
      break;
    }
  }
  else
  {
    if( MPL::rank() == 0 )
      std::cout << "skipping tests with 5 mpi tasks!" << std::endl;
  }

  actions::build_edges(*m);
  actions::build_median_dual_mesh(*m);

  BOOST_CHECK_CLOSE( test::dual_volume(*m), 2.*M_PI*M_PI, 1e-6 );

  std::stringstream filename; filename << "small_halo_p" << MPL::rank() << ".msh";
  Gmsh().write(*m,filename.str());
}
#endif

#if 1
BOOST_AUTO_TEST_CASE( test_t63 )
{
//  Mesh::Ptr m = test::generate_mesh( T63() );

  int nlat = 5;
  int lon[5] = {10, 12, 14, 16, 16};
  Mesh::Ptr m = test::generate_mesh(nlat, lon);

  actions::build_nodes_parallel_fields(m->function_space("nodes"));
  actions::build_periodic_boundaries(*m);
  actions::build_halo(*m,1);
  //actions::build_edges(*m);
  //actions::build_pole_edges(*m);
  //actions::build_edges_parallel_fields(m->function_space("edges"),m->function_space("nodes"));
  //actions::build_centroid_dual_mesh(*m);
  actions::renumber_nodes_glb_idx(m->function_space("nodes"));

  std::stringstream filename; filename << "T63_halo.msh";
  Gmsh().write(*m,filename.str());

//  BOOST_CHECK_CLOSE( test::dual_volume(*m), 2.*M_PI*M_PI, 1e-6 );

//  FunctionSpace& nodes = m->function_space("nodes");
//  FunctionSpace& edges = m->function_space("edges");
//  ArrayView<double,1> dual_volumes  ( nodes.field( "dual_volumes" ) );
//  ArrayView<double,2> dual_normals  ( edges.field( "dual_normals" ) );

//  std::string checksum;
//  checksum = nodes.checksum()->execute(dual_volumes);
//  DEBUG("dual_volumes checksum "<<checksum,0);

//  checksum = edges.checksum()->execute(dual_normals);
//  DEBUG("dual_normals checksum "<<checksum,0);

}
#endif

