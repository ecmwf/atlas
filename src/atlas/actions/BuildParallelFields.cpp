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
#include "atlas/mpl/MPL.hpp"
#include "atlas/mpl/GatherScatter.hpp"

namespace atlas {
namespace actions {

FieldT<int>& build_nodes_partition ( FunctionSpace& nodes );
FieldT<int>& build_nodes_remote_idx( FunctionSpace& nodes );
FieldT<int>& build_nodes_global_idx( FunctionSpace& nodes );
FieldT<int>& build_edges_partition ( FunctionSpace& edges, FunctionSpace& nodes );
FieldT<int>& build_edges_remote_idx( FunctionSpace& edges, FunctionSpace& nodes );
FieldT<int>& build_edges_global_idx( FunctionSpace& edges, FunctionSpace& nodes );

// ------------------------------------------------------------------

namespace {

struct Node
{
  Node(int gid, int idx)
  {
    g = gid;
    i = idx;
  }
  int g,i;
  bool operator < (const Node& other) const
  {
    return ( g<other.g );
  }
};

}

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

  build_nodes_partition ( nodes );
  build_nodes_remote_idx( nodes );
  build_nodes_global_idx( nodes );

}

// ------------------------------------------------------------------

void build_edges_parallel_fields( FunctionSpace& edges, FunctionSpace& nodes )
{
  ASSERT( nodes.has_field("coordinates") );
  ASSERT( nodes.has_field("partition") );
  ASSERT( nodes.has_field("remote_idx") );
  ASSERT( nodes.has_field("glb_idx") );

  build_edges_partition ( edges, nodes );
  build_edges_remote_idx( edges, nodes );
  build_edges_global_idx( edges, nodes );

  ASSERT( edges.has_field("partition") );
  ASSERT( edges.has_field("remote_idx") );
  ASSERT( edges.has_field("glb_idx") );

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

void renumber_nodes_glb_idx( FunctionSpace& nodes )
{
  int mypart = MPL::rank();
  int nparts = MPL::size();
  int root = 0;

  if( ! nodes.has_field("glb_idx") )
  {
    ArrayView<int,1> glb_idx ( nodes.create_field<int>("glb_idx",1) );
    glb_idx = -1;
  }

  ArrayView<double,2> latlon  ( nodes.field("coordinates") );
  ArrayView<int,   1> glb_idx ( nodes.field("glb_idx"    ) );

  /*
   * Sorting following gidx will define global order of
   * gathered fields. Special care needs to be taken for
   * pole edges, as their centroid might coincide with
   * other edges
   */
  int nb_nodes = glb_idx.extents()[0];
  for( int jnode=0; jnode<nb_nodes; ++jnode )
  {
    if( glb_idx(jnode) <= 0 )
      glb_idx(jnode) = LatLonPoint( latlon[jnode] ).uid();
  }

  /*
   * REMOTE INDEX BASE = 1
   */
  const int ridx_base = 1;

  // 1) Gather all global indices, together with location
  Array<int> loc_id_arr(nb_nodes);
  ArrayView<int,1> loc_id(loc_id_arr);

  for( int jnode=0; jnode<nb_nodes; ++jnode )
  {
    loc_id(jnode) = glb_idx(jnode);
  }

  std::vector<int> recvcounts(MPL::size());
  std::vector<int> recvdispls(MPL::size());
  MPL_CHECK_RESULT( MPI_Gather( &nb_nodes, 1, MPI_INT,
                                recvcounts.data(), 1, MPI_INT, root, MPI_COMM_WORLD) );
  recvdispls[0]=0;
  for (int jpart=1; jpart<nparts; ++jpart) // start at 1
  {
    recvdispls[jpart]=recvcounts[jpart-1]+recvdispls[jpart-1];
  }
  int glb_nb_nodes = std::accumulate(recvcounts.begin(),recvcounts.end(),0);

  Array<int> glb_id_arr(glb_nb_nodes);
  ArrayView<int,1> glb_id(glb_id_arr);

  MPL_CHECK_RESULT(
        MPI_Gatherv( loc_id.data(), nb_nodes, MPI_INT,
                     glb_id.data(), recvcounts.data(), recvdispls.data(), MPI_INT,
                     root, MPI_COMM_WORLD) );


  // 2) Sort all global indices, and renumber from 1 to glb_nb_edges
  std::vector<Node> node_sort; node_sort.reserve(glb_nb_nodes);
  for( int jnode=0; jnode<glb_id.extent(0); ++jnode )
  {
    node_sort.push_back( Node(glb_id(jnode),jnode) );
  }
  std::sort(node_sort.begin(), node_sort.end());

  // Assume edge gid start
  int gid=0;
  for( int jnode=0; jnode<node_sort.size(); ++jnode )
  {
    if( jnode == 0 )
    {
      ++gid;
    }
    else if( node_sort[jnode].g != node_sort[jnode-1].g )
    {
      ++gid;
    }
    int inode = node_sort[jnode].i;
    glb_id(inode) = gid;
  }

  // 3) Scatter renumbered back
  MPL_CHECK_RESULT(
        MPI_Scatterv( glb_id.data(), recvcounts.data(), recvdispls.data(), MPI_INT,
                      loc_id.data(), nb_nodes, MPI_INT,
                      root, MPI_COMM_WORLD) );

  for( int jnode=0; jnode<nb_nodes; ++jnode )
  {
    glb_idx(jnode) = loc_id(jnode);
  }
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


  int varsize=3;

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
        Extents(recv_needed[ proc[jpart] ].size()/varsize,varsize) );
    for( int jnode=0; jnode<recv_node.extents()[0]; ++jnode )
    {
      LatLonPoint ll( recv_node[jnode] );
      if( lookup.count(ll) )
      {
        send_found[ proc[jpart] ].push_back( recv_node(jnode,2) );
        send_found[ proc[jpart] ].push_back( lookup[ll] );
      }
      else
      {
        std::stringstream msg;
        msg << "[" << MPL::rank() << "] " << "Node requested by rank ["<<jpart<<"] with coords ";
        msg << "("<<ll.x*180./M_PI*1.e-6 << ","<<ll.y*180./M_PI*1.e-6<<") that should be owned is not found";
        throw eckit::SeriousBug(msg.str(),Here());
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

FieldT<int>& build_edges_partition( FunctionSpace& edges, FunctionSpace& nodes )
{
  int mypart = MPL::rank();
  int nparts = MPL::size();

  if( ! edges.has_field("partition") ) edges.create_field<int>("partition",1) ;
  ArrayView<int,1> edge_part  ( edges.field("partition") );
  IndexView<int,2> edge_nodes ( edges.field("nodes")     );
  ArrayView<int,1> is_pole_edge;
  bool has_pole_edges = false;
  if( edges.has_field("is_pole_edge") )
  {
    has_pole_edges = true;
    is_pole_edge = ArrayView<int,1>( edges.field("is_pole_edge") );
  }

  ArrayView<int,1> node_part  ( nodes.field("partition") );
  ArrayView<double,2> latlon     ( nodes.field("coordinates") );

  // Set all edges partition to the minimum partition number of nodes
  int nb_edges = edges.extents()[0];
  for( int jedge=0; jedge<nb_edges; ++jedge )
  {
    int p1 = node_part( edge_nodes(jedge,0) );
    int p2 = node_part( edge_nodes(jedge,1) );
    edge_part(jedge) = std::min( p1,p2 );
  }

  // Some edges might not exist on the assigned part. Check this, and change
  // the partition numbers for those edges
  std::map<int,int> lookup;
  int varsize=2;
  double centroid[2];
  std::vector< std::vector<int> > send_needed( MPL::size() );
  std::vector< std::vector<int> > recv_needed( MPL::size() );
  for( int jedge=0; jedge<nb_edges; ++jedge )
  {
    centroid[XX] = 0.5*(latlon( edge_nodes(jedge,0), XX ) + latlon( edge_nodes(jedge,1), XX ) );
    centroid[YY] = 0.5*(latlon( edge_nodes(jedge,0), YY ) + latlon( edge_nodes(jedge,1), YY ) );
    if( has_pole_edges && is_pole_edge(jedge) )
    {
      centroid[YY] = centroid[YY] > 0 ? M_PI_2 : -M_PI_2;
    }
    LatLonPoint ll(centroid);
    bool needed(false);
    if     ( ll.x <   BC::WEST) { needed = true; while(ll.x <  BC::WEST) { ll.x += BC::EAST; } }
    else if( ll.x >=  BC::EAST) { needed = true; while(ll.x >= BC::EAST) { ll.x -= BC::EAST; } }
    int uid = ll.uid();
    if( edge_part(jedge)==mypart && !needed )
    {
      lookup[ uid ] = jedge;
    }
    else
    {
      send_needed[ edge_part(jedge) ].push_back( uid   );
      send_needed[ edge_part(jedge) ].push_back( jedge );
    }
  }

  MPL::Alltoall( send_needed, recv_needed );

  std::vector< std::vector<int> > send_notfound( MPL::size() );
  std::vector< std::vector<int> > recv_notfound( MPL::size() );

  for( int jpart=0; jpart<nparts; ++jpart )
  {
    ArrayView<int,2> recv_edge( recv_needed[ jpart ].data(),
        Extents(recv_needed[ jpart ].size()/varsize,varsize) );
    for( int jedge=0; jedge<recv_edge.extents()[0]; ++jedge )
    {
      int uid      = recv_edge(jedge,0);
      int recv_idx = recv_edge(jedge,1);
      if( lookup.count(uid) == 0 )
        send_notfound[ jpart ].push_back( recv_idx );
    }
  }

  MPL::Alltoall( send_notfound, recv_notfound );

  for( int jpart=0; jpart<nparts; ++jpart )
  {
    for( int jedge=0; jedge<recv_notfound[jpart].size(); ++jedge )
    {
      int iedge = recv_notfound[jpart][jedge];
      int p1 = node_part( edge_nodes(iedge,0) );
      int p2 = node_part( edge_nodes(iedge,1) );
      edge_part(iedge) = std::max(p1,p2);
      //DEBUG("change owner of " << edge_gidx(iedge) << " from " << std::min(p1,p2) << "to " << std::max(p1,p2));
    }
  }

  /// TODO: Make sure that the edge-partition is at least one of the partition numbers of the
  /// neighbouring elements.
  /// Because of this problem, the size of the halo should be set to 2 instead of 1!!!

  return edges.field<int>("partition");
}

FieldT<int>& build_edges_remote_idx( FunctionSpace& edges, FunctionSpace& nodes )
{
  int mypart = MPL::rank();
  int nparts = MPL::size();

  if( ! edges.has_field("remote_idx") ) ( edges.create_field<int>("remote_idx",1) );
  IndexView<int,   2> edge_nodes ( edges.field("nodes")       );
  IndexView<int,   1> edge_ridx  ( edges.field("remote_idx")  );
  ArrayView<int,   1> edge_part  ( edges.field("partition")   );
  ArrayView<double,2> latlon     ( nodes.field("coordinates") );
  ArrayView<int,1> is_pole_edge;
  bool has_pole_edges = false;
  if( edges.has_field("is_pole_edge") )
  {
    has_pole_edges = true;
    is_pole_edge = ArrayView<int,1>( edges.field("is_pole_edge") );
  }

  const int nb_nodes = nodes.extents()[0];
  const int nb_edges = edges.extents()[0];

  int varsize=2;
  double centroid[2];
  std::vector< std::vector<int> > send_needed( MPL::size() );
  std::vector< std::vector<int> > recv_needed( MPL::size() );
  int sendcnt=0;
  std::map<int,int> lookup;

  for( int jedge=0; jedge<nb_edges; ++jedge )
  {
    centroid[XX] = 0.5*(latlon( edge_nodes(jedge,0), XX ) + latlon( edge_nodes(jedge,1), XX ) );
    centroid[YY] = 0.5*(latlon( edge_nodes(jedge,0), YY ) + latlon( edge_nodes(jedge,1), YY ) );
    if( has_pole_edges && is_pole_edge(jedge) )
    {
      centroid[YY] = centroid[YY] > 0 ? M_PI_2 : -M_PI_2;
    }
    LatLonPoint ll(centroid);
    bool needed(false);
    if     ( ll.x <   BC::WEST) { needed = true; while(ll.x <  BC::WEST) { ll.x += BC::EAST; } }
    else if( ll.x >=  BC::EAST) { needed = true; while(ll.x >= BC::EAST) { ll.x -= BC::EAST; } }
    int uid = ll.uid();
    if( edge_part(jedge)==mypart && !needed )
    {
      lookup[ uid ] = jedge;
      edge_ridx(jedge) = jedge;
    }
    else
    {
      send_needed[ edge_part(jedge) ].push_back( uid   );
      send_needed[ edge_part(jedge) ].push_back( jedge );
      sendcnt++;
    }
  }

  MPL::Alltoall( send_needed, recv_needed );

  std::vector< std::vector<int> > send_found( MPL::size() );
  std::vector< std::vector<int> > recv_found( MPL::size() );

  std::map<int,int>::iterator found;
  for( int jpart=0; jpart<nparts; ++jpart )
  {
    ArrayView<int,2> recv_edge( recv_needed[ jpart ].data(),
        Extents(recv_needed[ jpart ].size()/varsize,varsize) );
    for( int jedge=0; jedge<recv_edge.extents()[0]; ++jedge )
    {
      int recv_uid = recv_edge(jedge,0);
      int recv_idx = recv_edge(jedge,1);
      found = lookup.find(recv_uid);
      if( found != lookup.end() )
      {
        send_found[ jpart ].push_back( recv_idx );
        send_found[ jpart ].push_back( found->second );
      }
      else
      {
        std::stringstream msg;
        msg << "[" << MPL::rank() << "] " << "Edge with uid " << recv_uid << " requested by rank ["<<jpart<<"]";
        msg << " that should be owned is not found";
        throw eckit::SeriousBug(msg.str(),Here());
      }
    }
  }

  MPL::Alltoall( send_found, recv_found );

  for( int jpart=0; jpart<nparts; ++jpart )
  {
    ArrayView<int,2> recv_edge( recv_found[ jpart ].data(),
        Extents(recv_found[ jpart ].size()/2,2) );
    for( int jedge=0; jedge<recv_edge.extents()[0]; ++jedge )
    {
      edge_ridx( recv_edge(jedge,0) ) = recv_edge(jedge,1);
    }
  }
  return edges.field<int>("remote_idx");

}

FieldT<int>& build_edges_global_idx( FunctionSpace& edges, FunctionSpace& nodes )
{
  int mypart = MPL::rank();
  int nparts = MPL::size();
  int root = 0;

  if( ! edges.has_field("glb_idx") )
  {
    ArrayView<int,1> edge_gidx ( edges.create_field<int>("glb_idx",1) );
    edge_gidx = -1;
  }

  ArrayView<int,   1> edge_gidx  ( edges.field("glb_idx")     );
  IndexView<int,   2> edge_nodes ( edges.field("nodes")       );
  ArrayView<double,2> latlon     ( nodes.field("coordinates") );
  ArrayView<int,1> is_pole_edge;
  bool has_pole_edges = false;
  if( edges.has_field("is_pole_edge") )
  {
    has_pole_edges = true;
    is_pole_edge = ArrayView<int,1>( edges.field("is_pole_edge") );
  }

  /*
   * Sorting following edge_gidx will define global order of
   * gathered fields. Special care needs to be taken for
   * pole edges, as their centroid might coincide with
   * other edges
   */
  double centroid[2];
  int nb_edges = edges.extents()[0];
  for( int jedge=0; jedge<nb_edges; ++jedge )
  {
    if( edge_gidx(jedge) <= 0 )
    {
      centroid[XX] = 0.5*(latlon( edge_nodes(jedge,0), XX ) + latlon( edge_nodes(jedge,1), XX ) );
      centroid[YY] = 0.5*(latlon( edge_nodes(jedge,0), YY ) + latlon( edge_nodes(jedge,1), YY ) );
      if( has_pole_edges && is_pole_edge(jedge) )
      {
        centroid[YY] = centroid[YY] > 0 ? M_PI_2 : -M_PI_2;
      }
      LatLonPoint ll(centroid);
      edge_gidx(jedge) = ll.uid();
    }
  }

  /*
   * REMOTE INDEX BASE = 1
   */
  const int ridx_base = 1;

  // 1) Gather all global indices, together with location
  Array<int> loc_edge_id_arr(nb_edges);
  ArrayView<int,1> loc_edge_id(loc_edge_id_arr);

  for( int jedge=0; jedge<nb_edges; ++jedge )
  {
    loc_edge_id(jedge) = edge_gidx(jedge);
  }

  std::vector<int> recvcounts(MPL::size());
  std::vector<int> recvdispls(MPL::size());
  MPL_CHECK_RESULT( MPI_Gather( &nb_edges, 1, MPI_INT,
                                recvcounts.data(), 1, MPI_INT, root, MPI_COMM_WORLD) );
  recvdispls[0]=0;
  for (int jpart=1; jpart<nparts; ++jpart) // start at 1
  {
    recvdispls[jpart]=recvcounts[jpart-1]+recvdispls[jpart-1];
  }
  int glb_nb_edges = std::accumulate(recvcounts.begin(),recvcounts.end(),0);

  Array<int> glb_edge_id_arr(glb_nb_edges);
  ArrayView<int,1> glb_edge_id(glb_edge_id_arr);

  MPL_CHECK_RESULT(
        MPI_Gatherv( loc_edge_id.data(), nb_edges, MPI_INT,
                     glb_edge_id.data(), recvcounts.data(), recvdispls.data(), MPI_INT,
                     root, MPI_COMM_WORLD) );


  // 2) Sort all global indices, and renumber from 1 to glb_nb_edges
  std::vector<Node> edge_sort; edge_sort.reserve(glb_nb_edges);
  for( int jedge=0; jedge<glb_edge_id.extent(0); ++jedge )
  {
    edge_sort.push_back( Node(glb_edge_id(jedge),jedge) );
  }
  std::sort(edge_sort.begin(), edge_sort.end());

  // Assume edge gid start
  int gid=200000;
  for( int jedge=0; jedge<edge_sort.size(); ++jedge )
  {
    if( jedge == 0 )
    {
      ++gid;
    }
    else if( edge_sort[jedge].g != edge_sort[jedge-1].g )
    {
      ++gid;
    }
    int iedge = edge_sort[jedge].i;
    glb_edge_id(iedge) = gid;
  }

  // 3) Scatter renumbered back
  MPL_CHECK_RESULT(
        MPI_Scatterv( glb_edge_id.data(), recvcounts.data(), recvdispls.data(), MPI_INT,
                      loc_edge_id.data(), nb_edges, MPI_INT,
                      root, MPI_COMM_WORLD) );



  for( int jedge=0; jedge<nb_edges; ++jedge )
  {
    edge_gidx(jedge) = loc_edge_id(jedge);
  }

  return edges.field<int>("glb_idx");
}

// ------------------------------------------------------------------
// C wrapper interfaces to C++ routines

void atlas__build_parallel_fields ( Mesh* mesh) {
  build_parallel_fields(*mesh);
}
void atlas__build_nodes_parallel_fields (FunctionSpace* nodes) {
  build_nodes_parallel_fields(*nodes);
}
void atlas__build_edges_parallel_fields (FunctionSpace* edges, FunctionSpace* nodes) {
  build_edges_parallel_fields(*edges, *nodes);
}
void atlas__renumber_nodes_glb_idx (FunctionSpace* nodes)
{
  renumber_nodes_glb_idx(*nodes);
}


// ------------------------------------------------------------------



} // namespace actions
} // namespace atlas

