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

#define BOOST_TEST_MODULE TestDistributeMesh
#define BOOST_UNIT_TEST_FRAMEWORK_HEADER_ONLY
#include "ecbuild/boost_test_framework.h"

#include "atlas/mpl/MPL.hpp"
#include "atlas/atlas_config.h"
#include "tests/TestMeshes.hpp"
#include "atlas/io/Gmsh.hpp"
#include "atlas/mesh/Mesh.hpp"
#include "atlas/mesh/FunctionSpace.hpp"
#include "atlas/util/IndexView.hpp"
#include "atlas/util/ArrayView.hpp"
#include "atlas/actions/BuildHalo.hpp"
#include "atlas/actions/BuildParallelFields.hpp"
#include "atlas/actions/BuildPeriodicBoundaries.hpp"
#include "atlas/actions/BuildEdges.hpp"
#include "atlas/actions/BuildDualMesh.hpp"
#include "atlas/actions/DistributeMesh.hpp"
#include "atlas/actions/WriteLoadBalanceReport.hpp"
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

BOOST_AUTO_TEST_CASE( test_distribute_t63 )
{
  // Every task builds full mesh
  meshgen::RGGMeshGenerator generator;
  generator.options.set("nb_parts",1);
  generator.options.set("part",0);

      // int lon[] = {4,6,8,8,8};
      // test::TestGrid grid(5,lon);

      //  GG grid(120,60);
  T63 grid;


  generator.options.set("nb_parts", MPL::size());
  generator.options.set("part", MPL::rank());

  Mesh::Ptr m( generator.generate( grid ) );

      //  Mesh::Ptr m( Gmsh::read("unstr.msh") );

//  actions::distribute_mesh(*m);

  actions::build_parallel_fields(*m);
  actions::build_periodic_boundaries(*m);
  actions::build_halo(*m,2);
  actions::build_edges(*m);
  actions::build_pole_edges(*m);
  actions::build_median_dual_mesh(*m);
  BOOST_CHECK_CLOSE( test::dual_volume(*m), 2.*M_PI*M_PI, 0.0001 );
  double difference = 2.*M_PI*M_PI - test::dual_volume(*m);
  if( difference > 1e-8 )
  {
    std::cout << "difference = " << difference << std::endl;
  }

  std::stringstream filename; filename << "T63_dist_p" << MPL::rank() << ".msh";
  Gmsh().write(*m,filename.str());

  actions::write_load_balance_report(*m,"load_balance.dat");
}


