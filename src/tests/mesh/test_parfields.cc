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

#define BOOST_TEST_MODULE TestBuildParallelFields
#include "ecbuild/boost_test_framework.h"

#include "atlas/parallel/mpi/mpi.h"

#include "atlas/mesh/actions/BuildParallelFields.h"
#include "atlas/mesh/actions/BuildPeriodicBoundaries.h"
#include "atlas/internals/atlas_config.h"
#include "atlas/output/Gmsh.h"
#include "atlas/field/Field.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/util/Metadata.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/internals/Parameters.h"
#include "atlas/grid/partitioners/EqualRegionsPartitioner.h"
#include "atlas/grid/grids.h"
#include "atlas/mesh/generators/Structured.h"
#include "atlas/array/Array.h"
#include "atlas/array/ArrayView.h"
#include "atlas/array/IndexView.h"
#include "atlas/internals/Bitflags.h"

#include "tests/AtlasFixture.h"

using atlas::internals::Topology;

using namespace atlas::array;
using namespace atlas::internals;
using namespace atlas::output;
using namespace atlas::mesh::generators;

namespace atlas {
namespace test {

class IsGhost
{
public:
  IsGhost( const mesh::Nodes& nodes )
  {
    part_   = array::ArrayView<int,1> (nodes.partition() );
    ridx_   = IndexView<int,1> (nodes.remote_index() );
    mypart_ = parallel::mpi::comm().rank();
  }

  bool operator()(size_t idx) const
  {
    if( part_[idx] != mypart_ ) return true;
    if( ridx_[idx] != (int)idx     ) return true;
    return false;
  }
private:
  int mypart_;
  array::ArrayView<int,1> part_;
  IndexView<int,1> ridx_;
};


#define DISABLE if(0)
#define ENABLE if(1)

BOOST_GLOBAL_FIXTURE( AtlasFixture );

BOOST_AUTO_TEST_CASE( test1 )
{
  mesh::Mesh::Ptr m ( mesh::Mesh::create() );

  mesh::Nodes& nodes = m->nodes();
  nodes.resize(10);
  array::ArrayView<double,2> lonlat ( nodes.lonlat());
  array::ArrayView<gidx_t,1> glb_idx( nodes.global_index());
  array::ArrayView<int,1> part ( nodes.partition() );
  array::ArrayView<int,1> flags ( nodes.field("flags") );
  flags = Topology::NONE;

  // This is typically available
  glb_idx(0) = 1;    part(0) = 0;
  glb_idx(1) = 2;    part(1) = 0;
  glb_idx(2) = 3;    part(2) = 0;
  glb_idx(3) = 4;    part(3) = 0;
  glb_idx(4) = 5;    part(4) = 0;
  glb_idx(5) = 6;    part(5) = std::min(1,(int)parallel::mpi::comm().size()-1);
  glb_idx(6) = 7;    part(6) = std::min(1,(int)parallel::mpi::comm().size()-1);
  glb_idx(7) = 8;    part(7) = std::min(1,(int)parallel::mpi::comm().size()-1);
  glb_idx(8) = 9;    part(8) = std::min(1,(int)parallel::mpi::comm().size()-1);
  glb_idx(9) = 10;   part(9) = std::min(1,(int)parallel::mpi::comm().size()-1);

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

  mesh::actions::build_parallel_fields(*m);

  BOOST_REQUIRE( nodes.has_field("remote_idx") );

  IndexView<int,1> loc( nodes.remote_index() );
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

  test::IsGhost is_ghost( m->nodes() );

  switch ( parallel::mpi::comm().rank() )
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

  mesh::actions::build_periodic_boundaries(*m);

  // points 8 and 9 are periodic slave points of points 0 and 1
  BOOST_CHECK_EQUAL( part(8), 0 );
  BOOST_CHECK_EQUAL( part(9), 0 );
  BOOST_CHECK_EQUAL( loc(8), 0 );
  BOOST_CHECK_EQUAL( loc(9), 1 );
  if( parallel::mpi::comm().rank() == 1 )
  {
    BOOST_CHECK_EQUAL( is_ghost(8), true );
    BOOST_CHECK_EQUAL( is_ghost(9), true );
  }

}

BOOST_AUTO_TEST_CASE( test2 )
{
  util::Config meshgen_options;
  meshgen_options.set("angle",27.5);
  meshgen_options.set("triangulate",false);
  mesh::generators::Structured generate(meshgen_options);
  mesh::Mesh* m = generate(
      grid::gaussian::ClassicGaussian(32) );
  mesh::actions::build_parallel_fields(*m);

  mesh::Nodes& nodes = m->nodes();
  IndexView<int,1> loc_idx ( nodes.remote_index() );
  array::ArrayView<int,1> part    ( nodes.partition());
  array::ArrayView<gidx_t,1> glb_idx ( nodes.global_index() );

  test::IsGhost is_ghost(nodes);

  size_t nb_ghost = 0;
  for( size_t jnode=0; jnode<nodes.size(); ++jnode )
  {
    if( is_ghost(jnode) ) ++nb_ghost;
  }

  if( parallel::mpi::comm().rank() == 0 ) BOOST_CHECK_EQUAL( nb_ghost, 129 );
  if( parallel::mpi::comm().rank() == 1 ) BOOST_CHECK_EQUAL( nb_ghost, 0   );

  mesh::actions::build_periodic_boundaries(*m);

  int nb_periodic = -nb_ghost;
  for( size_t jnode=0; jnode<nodes.size(); ++jnode )
  {
    if( is_ghost(jnode) ) ++nb_periodic;
  }

  if( parallel::mpi::comm().rank() == 0 ) BOOST_CHECK_EQUAL( nb_periodic, 32 );
  if( parallel::mpi::comm().rank() == 1 ) BOOST_CHECK_EQUAL( nb_periodic, 32 );

  Gmsh("periodic.msh").write(*m);
  delete m;
}

} // namespace test
} // namespace atlas
