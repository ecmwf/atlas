/*
 * (C) Copyright 1996-2014 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */


#include <iostream>
#include <stdexcept>
#include <cmath>

#include "eckit/exception/Exceptions.h"
#include "atlas/mesh/Mesh.hpp"
#include "atlas/mesh/FunctionSpace.hpp"
#include "atlas/mesh/Field.hpp"
#include "atlas/actions/BuildParallelFields.hpp"
#include "atlas/mesh/Parameters.hpp"
#include "atlas/util/ArrayView.hpp"
#include "atlas/util/IndexView.hpp"
#include "atlas/util/Array.hpp"
#include "atlas/mesh/Util.hpp"

namespace atlas {
namespace actions {

FieldT<int>& build_nodes_global_idx( FunctionSpace& nodes );
FieldT<int>& build_nodes_partition ( FunctionSpace& nodes );
FieldT<int>& build_nodes_remote_idx( FunctionSpace& nodes );

// ------------------------------------------------------------------

void build_parallel_fields( Mesh& mesh )
{
  ASSERT( mesh.has_function_space("nodes") );

  build_nodes_parallel_fields( mesh.function_space("nodes") );
}

// ------------------------------------------------------------------

void build_nodes_parallel_fields( FunctionSpace& nodes )
{
  ASSERT( nodes.has_field("coordinates") );

  build_nodes_global_idx( nodes );
  build_nodes_partition ( nodes );
  build_nodes_remote_idx( nodes );
}

// ------------------------------------------------------------------

FieldT<int>& build_nodes_global_idx( FunctionSpace& nodes )
{
  if( ! nodes.has_field("glb_idx") )
  {
    ArrayView<int,1> glb_idx ( nodes.create_field<int>("glb_idx",1) );
    glb_idx = -1;
  }

  ArrayView<double,2> latlon  ( nodes.field("coordinates") );
  ArrayView<int,   1> glb_idx ( nodes.field("glb_idx"    ) );

  for( int jnode=0; jnode<glb_idx.extents()[0]; ++jnode )
  {
    if( glb_idx(jnode) <= 0 )
      glb_idx(jnode) = LatLonPoint( latlon[jnode] ).uid();
  }
  return nodes.field<int>("glb_idx");
}

// ------------------------------------------------------------------

FieldT<int>& build_nodes_remote_idx( FunctionSpace& nodes )
{
  int mypart = MPL::rank();
  int nparts = MPL::size();

  // This piece should be somewhere central ... could be NPROMA ?
  // ---------->
  std::vector< int > proc( nparts );
  for( int jpart=0; jpart<nparts; ++jpart )
    proc[jpart] = jpart;
  // <---------

  if( ! nodes.has_field("remote_idx") ) ( nodes.create_field<int>("remote_idx",1) );
  IndexView<int,   1> ridx   ( nodes.field("remote_idx")  );
  ArrayView<int,   1> part   ( nodes.field("partition")   );
  ArrayView<double,2> latlon ( nodes.field("coordinates") );
  int nb_nodes = nodes.extents()[0];

  std::vector< std::vector<int> > send_needed( MPL::size() );
  std::vector< std::vector<int> > recv_needed( MPL::size() );
  int sendcnt=0;
  std::map<LatLonPoint,int> lookup;
  for( int jnode=0; jnode<nb_nodes; ++jnode )
  {
    LatLonPoint ll(latlon[jnode]);
    if( part(jnode)==mypart )
    {
      lookup[ ll ] = jnode;
      ridx(jnode) = jnode;
    }
    else
    {
      send_needed[ proc[part(jnode)] ].push_back( ll.x  );
      send_needed[ proc[part(jnode)] ].push_back( ll.y  );
      send_needed[ proc[part(jnode)] ].push_back( jnode );
      sendcnt++;
    }
  }

  MPL::Alltoall( send_needed, recv_needed );

  std::vector< std::vector<int> > send_found( MPL::size() );
  std::vector< std::vector<int> > recv_found( MPL::size() );

  for( int jpart=0; jpart<nparts; ++jpart )
  {
    ArrayView<int,2> recv_node( recv_needed[ proc[jpart] ].data(),
        Extents(recv_needed[ proc[jpart] ].size()/3,3) );
    for( int jnode=0; jnode<recv_node.extents()[0]; ++jnode )
    {
      LatLonPoint ll( recv_node[jnode] );
      if( lookup.count(ll) )
      {
        send_found[ proc[jpart] ].push_back( recv_node(jnode,2) );
        send_found[ proc[jpart] ].push_back( lookup[ll] );
      }
    }
  }

  MPL::Alltoall( send_found, recv_found );

  for( int jpart=0; jpart<nparts; ++jpart )
  {
    ArrayView<int,2> recv_node( recv_found[ proc[jpart] ].data(),
        Extents(recv_found[ proc[jpart] ].size()/2,2) );
    for( int jnode=0; jnode<recv_node.extents()[0]; ++jnode )
    {
      ridx( recv_node(jnode,0) ) = recv_node(jnode,1);
    }
  }
  return nodes.field<int>("remote_idx");
}

// ------------------------------------------------------------------

FieldT<int>& build_nodes_partition( FunctionSpace& nodes )
{
  if( ! nodes.has_field("partition") )
  {
    ArrayView<int,1> part ( nodes.create_field<int>("partition",1) );
    part = MPL::rank();
  }
  return nodes.field<int>("partition");
}

// ------------------------------------------------------------------


// ------------------------------------------------------------------
// C wrapper interfaces to C++ routines

void atlas__build_parallel_fields ( Mesh* mesh) {
  build_parallel_fields(*mesh);
}
// ------------------------------------------------------------------



} // namespace actions
} // namespace atlas

