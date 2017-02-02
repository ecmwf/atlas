/*
 * (C) Copyright 1996-2017 ECMWF.
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
#include "ecbuild/boost_test_framework.h"

#include "atlas/parallel/mpi/mpi.h"
#include "atlas/atlas.h"
#include "tests/TestMeshes.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/output/Gmsh.h"
#include "atlas/array/ArrayView.h"
#include "atlas/array/MakeView.h"
#include "atlas/mesh/actions/BuildHalo.h"
#include "atlas/mesh/actions/BuildParallelFields.h"
#include "atlas/mesh/actions/BuildPeriodicBoundaries.h"
#include "atlas/mesh/actions/BuildEdges.h"
#include "atlas/mesh/actions/BuildDualMesh.h"
#include "atlas/mesh/actions/WriteLoadBalanceReport.h"
#include "atlas/internals/Parameters.h"
#include "atlas/grid/gaussian/classic/N.h"
#include "atlas/internals/IsGhost.h"
#include "atlas/runtime/Log.h"

#include "tests/AtlasFixture.h"


using namespace atlas;
using namespace atlas::output;
using namespace atlas::mesh::generators;
using namespace atlas::util;

namespace atlas {
namespace test {

double dual_volume(mesh::Mesh& mesh)
{
  mesh::Nodes& nodes = mesh.nodes();
  internals::IsGhost is_ghost_node(nodes);
  int nb_nodes = nodes.size();
  array::ArrayView<double,1> dual_volumes = array::make_view<double,1>( nodes.field("dual_volumes") );
  double area=0;

  for( int node=0; node<nb_nodes; ++node )
  {
    if( ! is_ghost_node(node) )
    {
      area += dual_volumes(node);
    }
  }



  parallel::mpi::comm().allReduceInPlace(area, eckit::mpi::sum());

  return area;
}

BOOST_GLOBAL_FIXTURE( AtlasFixture );

BOOST_AUTO_TEST_CASE( test_distribute_t63 )
{
  // Every task builds full mesh
//  mesh::generators::Structured generate( util::Config
//      ("nb_parts",1)
//      ("part",0) );
  mesh::generators::Structured generate;

      // long lon[] = {4,6,8,8,8};
      // test::TestGrid grid(5,lon);

      //  GG grid(120,60);
  grid::gaussian::ClassicGaussian grid(16);


  mesh::Mesh::Ptr m( generate( grid ) );

  // mesh::Mesh::Id meshid = m->id();

//  actions::distribute_mesh(*m);

  mesh::actions::build_parallel_fields(*m);
  mesh::actions::build_periodic_boundaries(*m);
  mesh::actions::build_halo(*m,1);
  //mesh::actions::renumber_nodes_glb_idx(m->nodes());
  mesh::actions::build_edges(*m);
  mesh::actions::build_pole_edges(*m);
  mesh::actions::build_edges_parallel_fields(*m);
  mesh::actions::build_median_dual_mesh(*m);

  double computed_dual_volume = test::dual_volume(*m);
  BOOST_CHECK_CLOSE( computed_dual_volume, 360.*180., 0.0001 );
  double difference = 360.*180. - computed_dual_volume;
  if( difference > 1e-8 )
  {
    std::cout << "difference = " << difference << std::endl;
  }

  Gmsh("dist.msh").write(*m);

  mesh::actions::write_load_balance_report(*m,"load_balance.dat");

  // mesh::Mesh& mesh1 = mesh::Mesh::from_id(meshid);
  mesh::Mesh& mesh1 = *m;
  BOOST_CHECK( mesh1.nodes().size() == m->nodes().size() );

  const array::ArrayView<int,1> part  = array::make_view<int,1>( m->nodes().partition() );
  const array::ArrayView<int,1> ghost = array::make_view<int,1>( m->nodes().ghost() );
  const array::ArrayView<int,1> flags = array::make_view<int,1>( m->nodes().field("flags") );

  Log::info() << "partition = [ ";
  for( size_t jnode=0; jnode<part.size(); ++jnode )
  {
    Log::info() << part(jnode) << " ";
  }
  Log::info() << "]" << std::endl;

  Log::info() << "ghost     = [ ";
  for( size_t jnode=0; jnode<part.size(); ++jnode )
  {
    Log::info() << ghost(jnode) << " ";
  }
  Log::info() << "]" << std::endl;

  Log::info() << "flags     = [ ";
  for( size_t jnode=0; jnode<part.size(); ++jnode )
  {
    Log::info() << internals::Topology::check(flags(jnode),internals::Topology::GHOST) << " ";
    BOOST_CHECK_EQUAL( internals::Topology::check(flags(jnode),internals::Topology::GHOST), ghost(jnode) );
  }
  Log::info() << "]" << std::endl;

}

} // namespace test
} // namespace atlas

