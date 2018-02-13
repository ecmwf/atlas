/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <sstream>
#include <algorithm>

#include "atlas/parallel/mpi/mpi.h"

#include "atlas/mesh/actions/BuildParallelFields.h"
#include "atlas/mesh/actions/BuildPeriodicBoundaries.h"
#include "atlas/library/config.h"
#include "atlas/output/Gmsh.h"
#include "atlas/field/Field.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/util/Metadata.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/util/CoordinateEnums.h"
#include "atlas/grid/detail/partitioner/EqualRegionsPartitioner.h"
#include "atlas/grid.h"
#include "atlas/meshgenerator/StructuredMeshGenerator.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/array.h"

#include "tests/AtlasTestEnvironment.h"

using Topology = atlas::mesh::Nodes::Topology;

using namespace atlas::array;
using namespace atlas::output;
using namespace atlas::meshgenerator;

namespace atlas {
namespace test {

//-----------------------------------------------------------------------------

class IsGhost
{
public:
  IsGhost( const mesh::Nodes& nodes ) :
    part_( make_view<int,1>(nodes.partition()) ),
    ridx_( make_indexview<int,1>(nodes.remote_index()) ),
    mypart_(parallel::mpi::comm().rank())
  {

  }

  bool operator()(size_t idx) const
  {
    if( part_(idx) != mypart_ ) return true;
    if( ridx_(idx) != (int)idx     ) return true;
    return false;
  }
private:
  array::ArrayView<int,1> part_;
  IndexView<int,1> ridx_;
  int mypart_;
};


#define DISABLE if(0)
#define ENABLE if(1)

//-----------------------------------------------------------------------------

CASE( "test1" )
{
  Mesh m;

  mesh::Nodes& nodes = m.nodes();
  nodes.resize(10);
  array::ArrayView<double,2> xy       = make_view<double,2>( nodes.xy());
  array::ArrayView<gidx_t,1> glb_idx  = make_view<gidx_t,1>( nodes.global_index());
  array::ArrayView<int   ,1> part     = make_view<int   ,1>( nodes.partition() );
  array::ArrayView<int   ,1> flags    = make_view<int   ,1>( nodes.field("flags") );
  flags.assign(Topology::NONE);

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

  xy(0,XX) = 0.;    xy(0,YY) = 80.;    Topology::set( flags(0), Topology::BC|Topology::WEST );
  xy(1,XX) = 0.;    xy(1,YY) =-80.;    Topology::set( flags(1), Topology::BC|Topology::WEST );
  xy(2,XX) = 90.;   xy(2,YY) = 80.;
  xy(3,XX) = 90.;   xy(3,YY) =-80.;
  xy(4,XX) = 180.;  xy(4,YY) = 80.;
  xy(5,XX) = 180.;  xy(5,YY) =-80.;
  xy(6,XX) = 270.;  xy(6,YY) = 80.;
  xy(7,XX) = 270.;  xy(7,YY) =-80.;
  xy(8,XX) = 360.;  xy(8,YY) = 80.;    Topology::set( flags(8), Topology::BC|Topology::EAST );
  xy(9,XX) = 360.;  xy(9,YY) =-80.;    Topology::set( flags(9), Topology::BC|Topology::EAST );

  mesh::actions::build_parallel_fields(m);

  EXPECT( nodes.has_field("remote_idx") );

  IndexView<int,1> loc = make_indexview<int,1>( nodes.remote_index() );
  EXPECT( loc(0) == 0 );
  EXPECT( loc(1) == 1 );
  EXPECT( loc(2) == 2 );
  EXPECT( loc(3) == 3 );
  EXPECT( loc(4) == 4 );
  EXPECT( loc(5) == 5 );
  EXPECT( loc(6) == 6 );
  EXPECT( loc(7) == 7 );
  EXPECT( loc(8) == 8 );
  EXPECT( loc(9) == 9 );

  test::IsGhost is_ghost( m.nodes() );

  switch ( parallel::mpi::comm().rank() )
  {
  case 0:
    EXPECT( is_ghost(0) == false );
    EXPECT( is_ghost(1) == false );
    EXPECT( is_ghost(2) == false );
    EXPECT( is_ghost(3) == false );
    EXPECT( is_ghost(4) == false );
    EXPECT( is_ghost(5) == true );
    EXPECT( is_ghost(6) == true );
    EXPECT( is_ghost(7) == true );
    EXPECT( is_ghost(8) == true );
    EXPECT( is_ghost(9) == true );
    break;
  case 1:
    EXPECT( is_ghost(0) == true );
    EXPECT( is_ghost(1) == true );
    EXPECT( is_ghost(2) == true );
    EXPECT( is_ghost(3) == true );
    EXPECT( is_ghost(4) == true );
    EXPECT( is_ghost(5) == false );
    EXPECT( is_ghost(6) == false );
    EXPECT( is_ghost(7) == false );
    EXPECT( is_ghost(8) == false );
    EXPECT( is_ghost(9) == false );
    break;
  }

  mesh::actions::build_periodic_boundaries(m);

  // points 8 and 9 are periodic slave points of points 0 and 1
  EXPECT( part(8) == 0 );
  EXPECT( part(9) == 0 );
  EXPECT( loc(8) == 0 );
  EXPECT( loc(9) == 1 );
  if( parallel::mpi::comm().rank() == 1 )
  {
    EXPECT( is_ghost(8) == true );
    EXPECT( is_ghost(9) == true );
  }

}

CASE( "test2" )
{
  util::Config meshgen_options;
  meshgen_options.set("angle",27.5);
  meshgen_options.set("triangulate",false);
  meshgenerator::StructuredMeshGenerator generate(meshgen_options);
  Mesh m = generate( Grid("N32") );
  mesh::actions::build_parallel_fields(m);

  mesh::Nodes& nodes = m.nodes();

  test::IsGhost is_ghost(nodes);

  size_t nb_ghost = 0;
  for( size_t jnode=0; jnode<nodes.size(); ++jnode )
  {
    if( is_ghost(jnode) ) ++nb_ghost;
  }

  ATLAS_DEBUG_VAR( nb_ghost );
  if( parallel::mpi::comm().rank() == 0 ) EXPECT( nb_ghost == 128 ); // South boundary of Northern hemisphere
  if( parallel::mpi::comm().rank() == 1 ) EXPECT( nb_ghost == 0   ); // Southern hemisphere has no ghosts

  mesh::actions::build_periodic_boundaries(m);

  int nb_periodic = -nb_ghost;
  for( size_t jnode=0; jnode<nodes.size(); ++jnode )
  {
    if( is_ghost(jnode) ) ++nb_periodic;
  }

  ATLAS_DEBUG_VAR( nb_periodic );

  if( parallel::mpi::comm().rank() == 0 ) EXPECT( nb_periodic == 33 ); // Periodic East boundary of Northern hemisphere (plus one point south)
  if( parallel::mpi::comm().rank() == 1 ) EXPECT( nb_periodic == 32 ); // Periodic East boundary of Southern hemisphere

  Gmsh("periodic.msh",util::Config("info",true)).write(m);
}

//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas


int main(int argc, char **argv) {
    return atlas::test::run( argc, argv );
}

