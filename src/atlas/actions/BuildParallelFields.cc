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
#include "eckit/log/Log.h"
#include "atlas/Mesh.h"
#include "atlas/FunctionSpace.h"
#include "atlas/Field.h"
#include "atlas/actions/BuildParallelFields.h"
#include "atlas/Parameters.h"
#include "atlas/util/ArrayView.h"
#include "atlas/util/IndexView.h"
#include "atlas/util/Array.h"
#include "atlas/Util.h"
#include "atlas/mpi/mpi.h"
#include "atlas/mpl/GatherScatter.h"

//#define DEBUGGING_PARFIELDS
#ifdef DEBUGGING_PARFIELDS
#define EDGE(jedge) "Edge("<<gidx(edge_nodes(jedge,0))<<"[p"<<node_part(edge_nodes(jedge,0))<<"] "<<gidx(edge_nodes(jedge,1))<<"p["<<node_part(edge_nodes(jedge,1))<<"])"
#define own1 2419089
#define own2 2423185
#define OWNED_EDGE(jedge) ((gidx(edge_nodes(jedge,0)) == own1 && gidx(edge_nodes(jedge,1)) == own2)\
                        || (gidx(edge_nodes(jedge,0)) == own2 && gidx(edge_nodes(jedge,1)) == own1))
#define per1 -1
#define per2 -1
#define PERIODIC_EDGE(jedge) ((gidx(edge_nodes(jedge,0)) == per1 && gidx(edge_nodes(jedge,1)) == per2)\
                          ||  (gidx(edge_nodes(jedge,0)) == per2 && gidx(edge_nodes(jedge,1)) == per1))
#define find1 -12
#define find2 -17
#define FIND_EDGE(jedge) ((gidx(edge_nodes(jedge,0)) == find1 && gidx(edge_nodes(jedge,1)) == find2)\
                      ||  (gidx(edge_nodes(jedge,0)) == find2 && gidx(edge_nodes(jedge,1)) == find1))
#define ownuid 547124520
#define OWNED_UID(UID) (UID == ownuid)
#endif

using eckit::Log;

namespace atlas {
namespace actions {

FieldT<int>& build_nodes_partition ( FunctionSpace& nodes );
FieldT<int>& build_nodes_remote_idx( FunctionSpace& nodes );
FieldT<gidx_t>& build_nodes_global_idx( FunctionSpace& nodes );
FieldT<int>& build_edges_partition ( FunctionSpace& edges, FunctionSpace& nodes );
FieldT<int>& build_edges_remote_idx( FunctionSpace& edges, FunctionSpace& nodes );
FieldT<gidx_t>& build_edges_global_idx( FunctionSpace& edges, FunctionSpace& nodes );

// ------------------------------------------------------------------

typedef gidx_t uid_t;

namespace {

struct Node
{
  Node(gidx_t gid, int idx)
  {
    g = gid;
    i = idx;
  }
  gidx_t g;
  gidx_t i;
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

FieldT<gidx_t>& build_nodes_global_idx( FunctionSpace& nodes )
{
  if( ! nodes.has_field("glb_idx") )
  {
    ArrayView<gidx_t,1> glb_idx ( nodes.create_field<gidx_t>("glb_idx",1) );
    glb_idx = -1;
  }

  ArrayView<gidx_t,1> glb_idx ( nodes.field("glb_idx"    ) );

  ComputeUid compute_uid(nodes);

  for( int jnode=0; jnode<glb_idx.shape(0); ++jnode )
  {
    if( glb_idx(jnode) <= 0 )
      glb_idx(jnode) = compute_uid(jnode);
  }
  return nodes.field<gidx_t>("glb_idx");
}

void renumber_nodes_glb_idx( FunctionSpace& nodes )
{
// TODO: ATLAS-14: fix renumbering of EAST periodic boundary points
// --> Those specific periodic points at the EAST boundary are not checked for uid,
//     and could receive different gidx for different tasks

  ComputeUid compute_uid(nodes);

  int mypart = mpi::rank();
  int nparts = mpi::size();
  int root = 0;

  if( ! nodes.has_field("glb_idx") )
  {
    ArrayView<gidx_t,1> glb_idx ( nodes.create_field<gidx_t>("glb_idx",1) );
    glb_idx = -1;
  }

  ArrayView<gidx_t,1> glb_idx ( nodes.field("glb_idx"    ) );

  /*
   * Sorting following gidx will define global order of
   * gathered fields. Special care needs to be taken for
   * pole edges, as their centroid might coincide with
   * other edges
   */
  int nb_nodes = glb_idx.shape(0);
  for( int jnode=0; jnode<nb_nodes; ++jnode )
  {
    if( glb_idx(jnode) <= 0 )
      glb_idx(jnode) = compute_uid(jnode);
  }


  // 1) Gather all global indices, together with location
  Array<uid_t> loc_id_arr(nb_nodes);
  ArrayView<uid_t,1> loc_id(loc_id_arr);

  for( int jnode=0; jnode<nb_nodes; ++jnode )
  {
    loc_id(jnode) = glb_idx(jnode);
  }

  std::vector<int> recvcounts(mpi::size());
  std::vector<int> recvdispls(mpi::size());
  ATLAS_MPI_CHECK_RESULT( MPI_Gather( &nb_nodes, 1, MPI_INT,
                                recvcounts.data(), 1, MPI_INT, root, mpi::Comm::instance()) );
  recvdispls[0]=0;
  for (int jpart=1; jpart<nparts; ++jpart) // start at 1
  {
    recvdispls[jpart]=recvcounts[jpart-1]+recvdispls[jpart-1];
  }
  int glb_nb_nodes = std::accumulate(recvcounts.begin(),recvcounts.end(),0);

  Array<uid_t> glb_id_arr(glb_nb_nodes);
  ArrayView<uid_t,1> glb_id(glb_id_arr);

  ATLAS_MPI_CHECK_RESULT(
        MPI_Gatherv( loc_id.data(), nb_nodes, mpi::datatype<uid_t>(),
                     glb_id.data(), recvcounts.data(), recvdispls.data(), mpi::datatype<uid_t>(),
                     root, mpi::Comm::instance()) );


  // 2) Sort all global indices, and renumber from 1 to glb_nb_edges
  std::vector<Node> node_sort; node_sort.reserve(glb_nb_nodes);
  for( int jnode=0; jnode<glb_id.shape(0); ++jnode )
  {
    node_sort.push_back( Node(glb_id(jnode),jnode) );
  }
  std::sort(node_sort.begin(), node_sort.end());

  // Assume edge gid start
  uid_t gid=0;
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
  ATLAS_MPI_CHECK_RESULT(
        MPI_Scatterv( glb_id.data(), recvcounts.data(), recvdispls.data(), mpi::datatype<uid_t>(),
                      loc_id.data(), nb_nodes, mpi::datatype<uid_t>(),
                      root, mpi::Comm::instance()) );

  for( int jnode=0; jnode<nb_nodes; ++jnode )
  {
    glb_idx(jnode) = loc_id(jnode);
  }
}

// ------------------------------------------------------------------

FieldT<int>& build_nodes_remote_idx( FunctionSpace& nodes )
{
  int mypart = mpi::rank();
  int nparts = mpi::size();

  ComputeUid compute_uid(nodes);

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
  int nb_nodes = nodes.shape(0);


  int varsize=2;

  std::vector< std::vector<uid_t> > send_needed( mpi::size() );
  std::vector< std::vector<uid_t> > recv_needed( mpi::size() );
  int sendcnt=0;
  std::map<uid_t,int> lookup;
  for( int jnode=0; jnode<nb_nodes; ++jnode )
  {
    uid_t uid = compute_uid(jnode);
    if( part(jnode)==mypart )
    {
      lookup[ uid ] = jnode;
      ridx(jnode) = jnode;
    }
    else
    {
      send_needed[ proc[part(jnode)] ].push_back( uid  );
      send_needed[ proc[part(jnode)] ].push_back( jnode );
      sendcnt++;
    }
  }

  mpi::all_to_all( send_needed, recv_needed );

  std::vector< std::vector<int> > send_found( mpi::size() );
  std::vector< std::vector<int> > recv_found( mpi::size() );

  for( int jpart=0; jpart<nparts; ++jpart )
  {
    ArrayView<uid_t,2> recv_node( recv_needed[ proc[jpart] ].data(),
        make_shape(recv_needed[ proc[jpart] ].size()/varsize,varsize) );
    for( int jnode=0; jnode<recv_node.shape(0); ++jnode )
    {
      uid_t uid = recv_node(jnode,0);
      int inode = recv_node(jnode,1);
      if( lookup.count(uid) )
      {
        send_found[ proc[jpart] ].push_back( inode );
        send_found[ proc[jpart] ].push_back( lookup[uid] );
      }
      else
      {
        std::stringstream msg;
        msg << "[" << mpi::rank() << "] " << "Node requested by rank ["<<jpart
            << "] with uid [" << uid << "] that should be owned is not found";
        throw eckit::SeriousBug(msg.str(),Here());
      }
    }
  }

  mpi::all_to_all( send_found, recv_found );

  for( int jpart=0; jpart<nparts; ++jpart )
  {
    ArrayView<int,2> recv_node( recv_found[ proc[jpart] ].data(),
        make_shape(recv_found[ proc[jpart] ].size()/2,2) );
    for( int jnode=0; jnode<recv_node.shape(0); ++jnode )
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
    part = mpi::rank();
  }
  return nodes.field<int>("partition");
}

// ------------------------------------------------------------------

FieldT<int>& build_edges_partition( FunctionSpace& edges, FunctionSpace& nodes )
{
  ComputeUid compute_uid(nodes);

  int mypart = mpi::rank();
  int nparts = mpi::size();

  if( ! edges.has_field("partition") ) edges.create_field<int>("partition",1) ;
  ArrayView<int,1> edge_part  ( edges.field("partition") );
  IndexView<int,2> edge_nodes ( edges.field("nodes")     );
  IndexView<int,2> edge_to_elem ( edges.field("to_elem")     );

  ArrayView<int,1> is_pole_edge;
  bool has_pole_edges = false;
  if( edges.has_field("is_pole_edge") )
  {
    has_pole_edges = true;
    is_pole_edge = ArrayView<int,1>( edges.field("is_pole_edge") );
  }

  ArrayView<int,1> node_part  ( nodes.field("partition") );
  ArrayView<double,2> latlon  ( nodes.field("coordinates") );
  ArrayView<int,   1> flags   ( nodes.field("flags")       );
#ifdef DEBUGGING_PARFIELDS
  ArrayView<gidx_t,   1> gidx    ( nodes.field("glb_idx")       );
#endif
  std::vector< IndexView<int,2> > elem_nodes( edges.mesh().nb_function_spaces() );
  std::vector< ArrayView<int,1> > elem_part ( edges.mesh().nb_function_spaces() );
  std::vector< ArrayView<gidx_t,1> > elem_glb_idx ( edges.mesh().nb_function_spaces() );

  for( int func_space_idx=0; func_space_idx<edges.mesh().nb_function_spaces(); ++func_space_idx)
  {
    FunctionSpace& elements = edges.mesh().function_space(func_space_idx);
    if( elements.metadata().get<int>("type") == Entity::ELEMS )
    {
      elem_nodes  [func_space_idx] = IndexView<int,2>( elements.field("nodes") );
      elem_part   [func_space_idx] = ArrayView<int,1>( elements.field("partition") );
      elem_glb_idx[func_space_idx] = ArrayView<gidx_t,1>( elements.field("glb_idx") );
//      int nb_elems = elem_nodes[func_space_idx].shape(0);
//      int nb_nodes_per_elem = elem_nodes[func_space_idx].shape(1);
//      for (int elem=0; elem<nb_elems; ++elem)
//      {
//        for (int n=0; n<nb_nodes_per_elem; ++n)
//        {
//          int node = elem_nodes[func_space_idx](elem,n);
//          node_to_elem[node].push_back( ElementRef(elements.index(),elem) );
//        }
//      }
    }
  }


  PeriodicTransform transform;

  int nb_edges = edges.shape(0);

  Array<int> periodic(nb_edges);

  for( int jedge=0; jedge<nb_edges; ++jedge )
  {
    periodic[jedge] = 0;
    int ip1 = edge_nodes(jedge,0);
    int ip2 = edge_nodes(jedge,1);
    int pn1 = node_part( ip1 );
    int pn2 = node_part( ip2 );
    if( pn1 == pn2 )
      edge_part(jedge) = pn1;
    else
    {
      if( edge_to_elem(jedge,2) < 0 ) // This is a edge at partition boundary
      {
        if( (Topology::check(flags(ip1),Topology::PERIODIC) && !Topology::check(flags(ip1),Topology::BC|Topology::WEST) &&
             Topology::check(flags(ip2),Topology::PERIODIC) && !Topology::check(flags(ip2),Topology::BC|Topology::WEST)) ||
            (Topology::check(flags(ip1),Topology::PERIODIC) && !Topology::check(flags(ip1),Topology::BC|Topology::WEST) &&
             Topology::check(flags(ip2),Topology::BC|Topology::WEST)) ||
            (Topology::check(flags(ip1),Topology::BC|Topology::WEST) &&
             Topology::check(flags(ip2),Topology::PERIODIC) && !Topology::check(flags(ip2),Topology::BC|Topology::WEST)) )
        {
          edge_part(jedge) = -1;
          if( Topology::check(flags(ip1),Topology::EAST ) )
            periodic[jedge] = -1;
          else
            periodic[jedge] = +1;
        }
        else if( Topology::check(flags(ip1),Topology::BC) && Topology::check(flags(ip2),Topology::BC) )
        {
          int pe1 = elem_part[ edge_to_elem(jedge,0) ](edge_to_elem(jedge,1));
          int pmx = std::max(pn1,pn2);
          if( pe1 == pmx )
            edge_part(jedge) = pmx;
          else
            edge_part(jedge) = std::min(pn1,pn2);
        }
        else
        {
          edge_part(jedge) = -1;
        }
      }
      else
      {
        int pe1 = elem_part[ edge_to_elem(jedge,0) ](edge_to_elem(jedge,1));
        int pe2 = elem_part[ edge_to_elem(jedge,2) ](edge_to_elem(jedge,3));
        int pc[] = {0,0};
        if( pn1 == pe1 ) ++pc[0];
        if( pn1 == pe2 ) ++pc[0];
        if( pn2 == pe1 ) ++pc[1];
        if( pn2 == pe2 ) ++pc[1];
        if     ( pc[0] > pc[1] ) edge_part(jedge) = pn1;
        else if( pc[0] < pc[1] ) edge_part(jedge) = pn2;
        else edge_part(jedge) = std::min(pn1,pn2);

#ifdef DEBUGGING_PARFIELDS
        if( OWNED_EDGE(jedge) )
          DEBUG( EDGE(jedge) << " -->    pn1,pn2,pe1,pe2 " << pn1 << "," << pn2<< "," << pe1 << "," << pe2);
#endif

      }
    }

#ifdef DEBUGGING_PARFIELDS
        if( OWNED_EDGE(jedge) )
          DEBUG( EDGE(jedge) << " -->    pn1,pn2 " << pn1 << "," << pn2);
#endif

#ifdef DEBUGGING_PARFIELDS_disable
    if( periodic[jedge] )
      DEBUG(EDGE(jedge)<<" is periodic  " << periodic[jedge]);
#endif

#ifdef DEBUGGING_PARFIELDS_disable
    if( PERIODIC_EDGE(jedge) )
      DEBUG( EDGE(jedge) <<": edge_part="<<edge_part(jedge)<<"  periodic="<<periodic[jedge]);
#endif
 }

  // In the periodic halo's, the partition MAY be assigned wrongly
  // following scoped piece of code will fix this.
  {
    std::map<uid_t,int> lookup;
    int varsize=2;
    double centroid[2];
    std::vector< std::vector<uid_t> > send_unknown( mpi::size() );
    std::vector< std::vector<uid_t> > recv_unknown( mpi::size() );
    for( int jedge=0; jedge<nb_edges; ++jedge )
    {
      int ip1 = edge_nodes(jedge,0);
      int ip2 = edge_nodes(jedge,1);
      int pn1 = node_part( ip1 );
      int pn2 = node_part( ip2 );

      centroid[XX] = 0.5*(latlon( ip1, XX ) + latlon( ip2, XX ) );
      centroid[YY] = 0.5*(latlon( ip1, YY ) + latlon( ip2, YY ) );
      if( has_pole_edges && is_pole_edge(jedge) )
      {
        centroid[YY] = centroid[YY] > 0 ? 90. : -90.;
      }

      transform(centroid,periodic[jedge]);
      uid_t uid = compute_uid(centroid);

      if( edge_part(jedge)==mypart )
      {
        lookup[ uid ] = jedge;
      }
      else
      {
        if( edge_part(jedge)<0 )
        {
          ASSERT(pn1!=pn2);
          send_unknown[ pn1 ].push_back( uid   );
          send_unknown[ pn1 ].push_back( jedge );
          send_unknown[ pn2 ].push_back( uid   );
          send_unknown[ pn2 ].push_back( jedge );

#ifdef DEBUGGING_PARFIELDS_disable
          if( PERIODIC_EDGE(jedge) )
            DEBUG_VAR( uid );
#endif

        }
      }
#ifdef DEBUGGING_PARFIELDS
        if( OWNED_EDGE(jedge) )
          DEBUG( EDGE(jedge) << " --> " << uid << "   part " << edge_part(jedge));
#endif

#ifdef DEBUGGING_PARFIELDS
        if( OWNED_UID(uid) )
        {
          double x1,y1, x2,y2, xe,ye;
          x1 = latlon(ip1,XX);
          y1 = latlon(ip1,YY);
          x2 = latlon(ip2,XX);
          y2 = latlon(ip2,YY);
          xe = centroid[XX];
          ye = centroid[YY];
          DEBUG( uid << " --> " << EDGE(jedge) << "   x1,y1 - x2,y2 - xe,ye " << x1<<","<<y1
                 << " - " << x2<<","<<y2<< " - " << xe <<","<<ye<< "     part " << edge_part(jedge));
        }
#endif

    }

    mpi::all_to_all( send_unknown, recv_unknown );

    // So now we have identified all possible edges with wrong partition.
    // We still need to check if it is actually wrong. This can be achieved
    // just by looking if the edge in question is owned by the possibly wrongly
    // assigned edge partition. If the edge is not found on that partition,
    // then its other node must be the edge partition.

    std::vector< std::vector<int> > send_found( mpi::size() );
    std::vector< std::vector<int> > recv_found( mpi::size() );

    for( int jpart=0; jpart<nparts; ++jpart )
    {
      ArrayView<uid_t,2> recv_edge( recv_unknown[ jpart ].data(),
          make_shape(recv_unknown[ jpart ].size()/varsize,varsize) );
      for( int jedge=0; jedge<recv_edge.shape(0); ++jedge )
      {
        uid_t uid      = recv_edge(jedge,0);
        int    recv_idx = recv_edge(jedge,1);
        if( lookup.count(uid) )
        {
          send_found[ jpart ].push_back( recv_idx );
        }
      }
    }

    mpi::all_to_all( send_found, recv_found );

    for( int jpart=0; jpart<nparts; ++jpart )
    {
      for( int jedge=0; jedge<recv_found[jpart].size(); ++jedge )
      {
        int iedge = recv_found[jpart][jedge];
#ifdef DEBUGGING_PARFIELDS
        //DEBUG(EDGE(iedge)<< " found on part " << jpart);
        if( edge_part(iedge) != -1 )
        {
          std::stringstream msg;
          msg << "Edge ("<<gidx(edge_nodes(iedge,0))<<"[p"<<node_part(edge_nodes(iedge,0))<<"] "<<gidx(edge_nodes(iedge,1))
              << "[p"<<node_part(edge_nodes(iedge,1))<<"]) from part " << jpart << " was already found on part " << edge_part(iedge);
          throw eckit::Exception(msg.str(),Here());
        }
#endif
        ASSERT( edge_part(iedge) == -1 );
        edge_part(iedge) = jpart;
      }
    }
  }
  for( int jedge=0; jedge<nb_edges; ++jedge )
  {
    int ip1 = edge_nodes(jedge,0);
    int ip2 = edge_nodes(jedge,1);
    if( edge_part(jedge) == -1 )
    {
      std::stringstream msg;
#ifdef DEBUGGING_PARFIELDS
      msg << EDGE(jedge) << " was not found on part "<< node_part(ip1) << " or " << node_part(ip2);
#else
      msg << "Edge was not found on part "<< node_part(ip1) << " or " <<node_part(ip2);
#endif
      throw eckit::SeriousBug(msg.str(),Here());
    }

#ifdef DEBUGGING_PARFIELDS
        if( OWNED_EDGE(jedge) )
          DEBUG( EDGE(jedge) << " -->    part " << edge_part(jedge));
#endif


#ifdef DEBUGGING_PARFIELDS_disable
    if( PERIODIC_EDGE(jedge) )
      DEBUG_VAR( "           the part is " << edge_part(jedge) );
#endif
  }
  /// TODO: Make sure that the edge-partition is at least one of the partition numbers of the
  /// neighbouring elements.
  /// Because of this problem, the size of the halo should be set to 2 instead of 1!!!
  /// This will be addressed with JIRA issue  ATLAS-12

  return edges.field<int>("partition");
}

FieldT<int>& build_edges_remote_idx( FunctionSpace& edges, FunctionSpace& nodes )
{
  ComputeUid compute_uid(nodes);

  int mypart = mpi::rank();
  int nparts = mpi::size();

  if( ! edges.has_field("remote_idx") ) ( edges.create_field<int>("remote_idx",1) );
  IndexView<int,   2> edge_nodes ( edges.field("nodes")       );
  IndexView<int,   1> edge_ridx  ( edges.field("remote_idx")  );
  ArrayView<int,   1> edge_part  ( edges.field("partition")   );
  ArrayView<double,2> latlon     ( nodes.field("coordinates") );
  ArrayView<int,   1> flags      ( nodes.field("flags")       );
#ifdef DEBUGGING_PARFIELDS
  ArrayView<gidx_t,   1> gidx      ( nodes.field("glb_idx")       );
  ArrayView<int,   1> node_part      ( nodes.field("partition")       );
#endif

  ArrayView<int,1> is_pole_edge;
  bool has_pole_edges = false;
  if( edges.has_field("is_pole_edge") )
  {
    has_pole_edges = true;
    is_pole_edge = ArrayView<int,1>( edges.field("is_pole_edge") );
  }

  const int nb_nodes = nodes.shape(0);
  const int nb_edges = edges.shape(0);

  double centroid[2];
  std::vector< std::vector<uid_t> > send_needed( mpi::size() );
  std::vector< std::vector<uid_t> > recv_needed( mpi::size() );
  int sendcnt=0;
  std::map<uid_t,int> lookup;

  PeriodicTransform transform;

  for( int jedge=0; jedge<nb_edges; ++jedge )
  {
    int ip1 = edge_nodes(jedge,0);
    int ip2 = edge_nodes(jedge,1);
    centroid[XX] = 0.5*(latlon( ip1, XX ) + latlon( ip2, XX ) );
    centroid[YY] = 0.5*(latlon( ip1, YY ) + latlon( ip2, YY ) );
    if( has_pole_edges && is_pole_edge(jedge) )
    {
      centroid[YY] = centroid[YY] > 0 ? 90. : -90.;
    }

    bool needed(false);

    if( (Topology::check(flags(ip1),Topology::PERIODIC) && !Topology::check(flags(ip1),Topology::BC|Topology::WEST) &&
         Topology::check(flags(ip2),Topology::PERIODIC) && !Topology::check(flags(ip2),Topology::BC|Topology::WEST)) ||
        (Topology::check(flags(ip1),Topology::PERIODIC) && !Topology::check(flags(ip1),Topology::BC|Topology::WEST) &&
         Topology::check(flags(ip2),Topology::BC|Topology::WEST)) ||
        (Topology::check(flags(ip1),Topology::BC|Topology::WEST) &&
         Topology::check(flags(ip2),Topology::PERIODIC) && !Topology::check(flags(ip2),Topology::BC|Topology::WEST)) )
    {
      needed = true;
      if( Topology::check(flags(ip1),Topology::EAST ) )
        transform(centroid,-1);
      else
        transform(centroid,+1);
    }

    uid_t uid = compute_uid(centroid);
    if( edge_part(jedge)==mypart && !needed ) // All interior edges fall here
    {
      lookup[ uid ] = jedge;
      edge_ridx(jedge) = jedge;

#ifdef DEBUGGING_PARFIELDS
      if( FIND_EDGE(jedge) )
      {
        DEBUG( "Found "<<EDGE(jedge));
      }
#endif

    }
    else // All ghost edges PLUS the periodic edges identified edges above fall here
    {
      send_needed[ edge_part(jedge) ].push_back( uid   );
      send_needed[ edge_part(jedge) ].push_back( jedge );
#ifdef DEBUGGING_PARFIELDS
      send_needed[ edge_part(jedge) ].push_back( gidx(ip1) );
      send_needed[ edge_part(jedge) ].push_back( gidx(ip2) );
      send_needed[ edge_part(jedge) ].push_back( node_part(ip1) );
      send_needed[ edge_part(jedge) ].push_back( node_part(ip2) );
#endif
      sendcnt++;
    }
  }
  int varsize=2;
#ifdef DEBUGGING_PARFIELDS
  varsize=6;
#endif
  mpi::all_to_all( send_needed, recv_needed );

  std::vector< std::vector<int> > send_found( mpi::size() );
  std::vector< std::vector<int> > recv_found( mpi::size() );

  std::map<uid_t,int>::iterator found;
  for( int jpart=0; jpart<nparts; ++jpart )
  {
    ArrayView<uid_t,2> recv_edge( recv_needed[ jpart ].data(),
        make_shape(recv_needed[ jpart ].size()/varsize,varsize) );
    for( int jedge=0; jedge<recv_edge.shape(0); ++jedge )
    {
      uid_t recv_uid = recv_edge(jedge,0);
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
#ifdef DEBUGGING_PARFIELDS
        msg << "Edge("<<recv_edge(jedge,2)<<"[p"<<recv_edge(jedge,4)<<"] "<<recv_edge(jedge,3)<<"[p"<<recv_edge(jedge,5)<<"])";
#else
        msg << "Edge with uid " << recv_uid;
#endif
        msg << " requested by rank ["<<jpart<<"]";
        msg << " that should be owned is not found. This could be because no halo was built.";
        //throw eckit::SeriousBug(msg.str(),Here());
        eckit::Log::warning() << msg.str() << " @ " << Here() << std::endl;
      }
    }
  }

  mpi::all_to_all( send_found, recv_found );

  for( int jpart=0; jpart<nparts; ++jpart )
  {
    ArrayView<int,2> recv_edge( recv_found[ jpart ].data(),
        make_shape(recv_found[ jpart ].size()/2,2) );
    for( int jedge=0; jedge<recv_edge.shape(0); ++jedge )
    {
      edge_ridx( recv_edge(jedge,0) ) = recv_edge(jedge,1);
    }
  }

  return edges.field<int>("remote_idx");

}

FieldT<gidx_t>& build_edges_global_idx( FunctionSpace& edges, FunctionSpace& nodes )
{
  ComputeUid compute_uid(nodes);

  int mypart = mpi::rank();
  int nparts = mpi::size();
  int root = 0;

  if( ! edges.has_field("glb_idx") )
  {
    ArrayView<gidx_t,1> edge_gidx ( edges.create_field<gidx_t>("glb_idx",1) );
    edge_gidx = -1;
  }

  ArrayView<gidx_t,1> edge_gidx  ( edges.field("glb_idx")     );
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
  int nb_edges = edges.shape(0);
  for( int jedge=0; jedge<nb_edges; ++jedge )
  {
    if( edge_gidx(jedge) <= 0 )
    {
      centroid[XX] = 0.5*(latlon( edge_nodes(jedge,0), XX ) + latlon( edge_nodes(jedge,1), XX ) );
      centroid[YY] = 0.5*(latlon( edge_nodes(jedge,0), YY ) + latlon( edge_nodes(jedge,1), YY ) );
      if( has_pole_edges && is_pole_edge(jedge) )
      {
        centroid[YY] = centroid[YY] > 0 ? 90. : -90.;
      }
      edge_gidx(jedge) = compute_uid(centroid);
    }
  }

  /*
   * REMOTE INDEX BASE = 1
   */
  const int ridx_base = 1;

  // 1) Gather all global indices, together with location
  Array<uid_t> loc_edge_id_arr(nb_edges);
  ArrayView<uid_t,1> loc_edge_id(loc_edge_id_arr);

  for( int jedge=0; jedge<nb_edges; ++jedge )
  {
    loc_edge_id(jedge) = edge_gidx(jedge);
  }

  std::vector<int> recvcounts(mpi::size());
  std::vector<int> recvdispls(mpi::size());
  ATLAS_MPI_CHECK_RESULT( MPI_Gather( &nb_edges, 1, MPI_INT,
                                recvcounts.data(), 1, MPI_INT, root, mpi::Comm::instance()) );
  recvdispls[0]=0;
  for (int jpart=1; jpart<nparts; ++jpart) // start at 1
  {
    recvdispls[jpart]=recvcounts[jpart-1]+recvdispls[jpart-1];
  }
  int glb_nb_edges = std::accumulate(recvcounts.begin(),recvcounts.end(),0);

  Array<uid_t> glb_edge_id_arr(glb_nb_edges);
  ArrayView<uid_t,1> glb_edge_id(glb_edge_id_arr);

  ATLAS_MPI_CHECK_RESULT(
        MPI_Gatherv( loc_edge_id.data(), nb_edges, mpi::datatype<uid_t>(),
                     glb_edge_id.data(), recvcounts.data(), recvdispls.data(), mpi::datatype<uid_t>(),
                     root, mpi::Comm::instance()) );


  // 2) Sort all global indices, and renumber from 1 to glb_nb_edges
  std::vector<Node> edge_sort; edge_sort.reserve(glb_nb_edges);
  for( int jedge=0; jedge<glb_edge_id.shape(0); ++jedge )
  {
    edge_sort.push_back( Node(glb_edge_id(jedge),jedge) );
  }
  std::sort(edge_sort.begin(), edge_sort.end());

  // Assume edge gid start
  uid_t gid=200000;
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
  ATLAS_MPI_CHECK_RESULT(
        MPI_Scatterv( glb_edge_id.data(), recvcounts.data(), recvdispls.data(), mpi::datatype<uid_t>(),
                      loc_edge_id.data(), nb_edges, mpi::datatype<uid_t>(),
                      root, mpi::Comm::instance()) );



  for( int jedge=0; jedge<nb_edges; ++jedge )
  {
    edge_gidx(jedge) = loc_edge_id(jedge);
  }

  return edges.field<gidx_t>("glb_idx");
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

