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
#include "ecbuild/boost_test_framework.h"

#include "atlas/mpi/mpi.h"
#include "atlas/atlas_config.h"
#include "tests/TestMeshes.h"
#include "atlas/io/Gmsh.h"
#include "atlas/Mesh.h"
#include "atlas/Nodes.h"
#include "atlas/FunctionSpace.h"
#include "atlas/util/IndexView.h"
#include "atlas/util/ArrayView.h"
#include "atlas/util/Array.h"
#include "atlas/actions/BuildHalo.h"
#include "atlas/actions/BuildParallelFields.h"
#include "atlas/actions/BuildPeriodicBoundaries.h"
#include "atlas/actions/BuildEdges.h"
#include "atlas/actions/BuildDualMesh.h"
#include "atlas/Parameters.h"
#include "atlas/util/IsGhost.h"


using namespace atlas::io;
using namespace atlas::util;
using namespace atlas::meshgen;

namespace atlas {
namespace test {

double dual_volume(Mesh& mesh)
{
  Nodes& nodes = mesh.nodes();
  IsGhost is_ghost_node(nodes);
  int nb_nodes = nodes.shape(0);
  ArrayView<double,1> dual_volumes ( nodes.field("dual_volumes") );
  ArrayView<gidx_t,1> glb_idx ( nodes.field("glb_idx") );
  double area=0;
  for( int node=0; node<nb_nodes; ++node )
  {
    if( ! is_ghost_node(node) )
    {
      area += dual_volumes(node);
    }
  }
  ECKIT_MPI_CHECK_RESULT( MPI_Allreduce( MPI_IN_PLACE, &area, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD ) );
  return area;
}


struct MPIFixture {
    MPIFixture()  { eckit::mpi::init(); }
    ~MPIFixture() { eckit::mpi::finalize(); }
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


  if( eckit::mpi::size() == 5 )
  {
    IndexView<int,1> ridx ( m->nodes().field("remote_idx") );
    ArrayView<gidx_t,1> gidx ( m->nodes().field("glb_idx") );

    switch( eckit::mpi::rank() ) // with 5 tasks
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
    if( eckit::mpi::rank() == 0 )
      std::cout << "skipping tests with 5 mpi tasks!" << std::endl;
  }

  actions::build_edges(*m);
  actions::build_median_dual_mesh(*m);

  BOOST_CHECK_CLOSE( test::dual_volume(*m), 2.*M_PI*M_PI, 1e-6 );

  std::stringstream filename; filename << "small_halo_p" << eckit::mpi::rank() << ".msh";
  Gmsh().write(*m,filename.str());
}
#endif

#if 1
BOOST_AUTO_TEST_CASE( test_t63 )
{
//  Mesh::Ptr m = test::generate_mesh( T63() );

  int nlat = 5;
  long lon[5] = {10, 12, 14, 16, 16};
  Mesh::Ptr m = test::generate_mesh(nlat, lon);

  actions::build_nodes_parallel_fields(m->nodes());
  actions::build_periodic_boundaries(*m);
  actions::build_halo(*m,1);
  //actions::build_edges(*m);
  //actions::build_pole_edges(*m);
  //actions::build_edges_parallel_fields(m->function_space("edges"),m->nodes());
  //actions::build_centroid_dual_mesh(*m);
  actions::renumber_nodes_glb_idx(m->nodes());

  std::stringstream filename; filename << "T63_halo.msh";
  Gmsh().write(*m,filename.str());

//  BOOST_CHECK_CLOSE( test::dual_volume(*m), 2.*M_PI*M_PI, 1e-6 );

//  Nodes& nodes = m->nodes();
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

} // namespace test
} // namespace atlas

