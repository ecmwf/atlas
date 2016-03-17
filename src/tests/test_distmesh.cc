/*
 * (C) Copyright 1996-2016 ECMWF.
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
#include "atlas/util/io/Gmsh.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/array/IndexView.h"
#include "atlas/array/ArrayView.h"
#include "atlas/mesh/actions/BuildHalo.h"
#include "atlas/mesh/actions/BuildParallelFields.h"
#include "atlas/mesh/actions/BuildPeriodicBoundaries.h"
#include "atlas/mesh/actions/BuildEdges.h"
#include "atlas/mesh/actions/BuildDualMesh.h"
#include "atlas/mesh/actions/WriteLoadBalanceReport.h"
#include "atlas/internals/Parameters.h"
#include "atlas/grid/global/gaussian/classic/N.h"
#include "atlas/internals/IsGhost.h"
#include "atlas/runtime/Log.h"

using namespace atlas;
using namespace atlas::util::io;
using namespace atlas::mesh::generators;
using namespace atlas::util;

namespace atlas {
namespace test {

double dual_volume(mesh::Mesh& mesh)
{
  mesh::Nodes& nodes = mesh.nodes();
  internals::IsGhost is_ghost_node(nodes);
  int nb_nodes = nodes.size();
  array::ArrayView<double,1> dual_volumes ( nodes.field("dual_volumes") );
  array::ArrayView<gidx_t,1> glb_idx ( nodes.global_index() );
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
     MPIFixture()  { atlas_init(); }
    ~MPIFixture()  { atlas_finalize(); }
};

BOOST_GLOBAL_FIXTURE( MPIFixture );

BOOST_AUTO_TEST_CASE( test_distribute_t63 )
{
  // Every task builds full mesh
  mesh::generators::ReducedGridMeshGenerator generate;
  generate.options.set("nb_parts",1);
  generate.options.set("part",0);

      // long lon[] = {4,6,8,8,8};
      // test::TestGrid grid(5,lon);

      //  GG grid(120,60);
  grid::predefined::rgg::N16 grid;


  generate.options.set("nb_parts", eckit::mpi::size());
  generate.options.set("part", eckit::mpi::rank());

  mesh::Mesh::Ptr m( generate( grid ) );

  mesh::Mesh::Id meshid = m->id();

      //  mesh::Mesh::Ptr m( Gmsh::read("unstr.msh") );

//  actions::distribute_mesh(*m);

  mesh::actions::build_parallel_fields(*m);
  mesh::actions::build_periodic_boundaries(*m);
  mesh::actions::build_halo(*m,1);
  //mesh::actions::renumber_nodes_glb_idx(m->nodes());
  mesh::actions::build_edges(*m);
  mesh::actions::build_pole_edges(*m);
  mesh::actions::build_edges_parallel_fields(*m);
  mesh::actions::build_median_dual_mesh(*m);
  BOOST_CHECK_CLOSE( test::dual_volume(*m), 360.*180., 0.0001 );
  double difference = 360.*180. - test::dual_volume(*m);
  if( difference > 1e-8 )
  {
    std::cout << "difference = " << difference << std::endl;
  }

  std::stringstream filename; filename << "N32_dist.msh";
  Gmsh().write(*m,filename.str());

  mesh::actions::write_load_balance_report(*m,"load_balance.dat");

  mesh::Mesh& mesh1 = mesh::Mesh::from_id(meshid);
  BOOST_CHECK(mesh1.grid().same(grid));



  const array::ArrayView<int,1> part( m->nodes().partition() );
  const array::ArrayView<int,1> ghost( m->nodes().ghost() );
  const array::ArrayView<int,1> flags( m->nodes().field("flags") );

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

