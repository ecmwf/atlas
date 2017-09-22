/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <cmath>
#include <iostream>
#include <stdexcept>
#include "eckit/exception/Exceptions.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/mesh/HybridElements.h"
#include "atlas/mesh/actions/BuildParallelFields.h"
#include "atlas/field/Field.h"
#include "atlas/util/CoordinateEnums.h"
#include "atlas/util/Unique.h"
#include "atlas/mesh/detail/PeriodicTransform.h"
#include "atlas/array/ArrayView.h"
#include "atlas/array/IndexView.h"
#include "atlas/array.h"
#include "atlas/runtime/Log.h"
#include "atlas/runtime/ErrorHandling.h"
#include "atlas/runtime/Timer.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/parallel/GatherScatter.h"

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

using Topology = atlas::mesh::Nodes::Topology;
using atlas::util::UniqueLonLat;
using atlas::mesh::detail::PeriodicTransform;

namespace atlas {
namespace mesh {
namespace actions {

Field& build_nodes_partition ( mesh::Nodes& nodes );
Field& build_nodes_remote_idx( mesh::Nodes& nodes );
Field& build_nodes_global_idx( mesh::Nodes& nodes );
Field& build_edges_partition ( Mesh& mesh );
Field& build_edges_remote_idx( Mesh& mesh );
Field& build_edges_global_idx( Mesh& mesh );

//----------------------------------------------------------------------------------------------------------------------

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

//----------------------------------------------------------------------------------------------------------------------


void build_parallel_fields( Mesh& mesh )
{
  Timer scope_timer(__FUNCTION__);
  build_nodes_parallel_fields( mesh.nodes() );
}

//----------------------------------------------------------------------------------------------------------------------

void build_nodes_parallel_fields( mesh::Nodes& nodes )
{
  Timer scope_timer(__FUNCTION__);
  bool parallel = false;
  nodes.metadata().get("parallel",parallel);
  if( ! parallel )
  {
    build_nodes_partition ( nodes );
    build_nodes_remote_idx( nodes );
    build_nodes_global_idx( nodes );
  }
  nodes.metadata().set("parallel",true);
}

//----------------------------------------------------------------------------------------------------------------------

void build_edges_parallel_fields( Mesh& mesh )
{
  Timer scope_timer(__FUNCTION__);
  build_edges_partition ( mesh );
  build_edges_remote_idx( mesh );
  build_edges_global_idx( mesh );
}

//----------------------------------------------------------------------------------------------------------------------

Field& build_nodes_global_idx( mesh::Nodes& nodes )
{
  Timer scope_timer(__FUNCTION__);

  array::ArrayView<gidx_t,1> glb_idx = array::make_view<gidx_t,1>( nodes.global_index() );

  UniqueLonLat compute_uid(nodes);

  for( size_t jnode=0; jnode<glb_idx.shape(0); ++jnode )
  {
    if( glb_idx(jnode) <= 0 )
      glb_idx(jnode) = compute_uid(jnode);
  }
  return nodes.global_index();
}

void renumber_nodes_glb_idx( mesh::Nodes& nodes )
{
  Timer scope_timer(__FUNCTION__);

// TODO: ATLAS-14: fix renumbering of EAST periodic boundary points
// --> Those specific periodic points at the EAST boundary are not checked for uid,
//     and could receive different gidx for different tasks

  UniqueLonLat compute_uid(nodes);

  // unused // int mypart = parallel::mpi::comm().rank();
  int nparts = parallel::mpi::comm().size();
  size_t root = 0;

  array::ArrayView<gidx_t,1> glb_idx = array::make_view<gidx_t,1> ( nodes.global_index() );

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
  array::ArrayT<uid_t> loc_id_arr( nb_nodes );
  array::ArrayView<uid_t,1> loc_id = array::make_view<uid_t,1>(loc_id_arr);

  for( int jnode=0; jnode<nb_nodes; ++jnode )
  {
    loc_id(jnode) = glb_idx(jnode);
  }

  std::vector<int> recvcounts(parallel::mpi::comm().size());
  std::vector<int> recvdispls(parallel::mpi::comm().size());

  parallel::mpi::comm().gather(nb_nodes, recvcounts, root);

  recvdispls[0]=0;
  for (int jpart=1; jpart<nparts; ++jpart) // start at 1
  {
    recvdispls[jpart]=recvcounts[jpart-1]+recvdispls[jpart-1];
  }
  int glb_nb_nodes = std::accumulate(recvcounts.begin(),recvcounts.end(),0);

  array::ArrayT<uid_t> glb_id_arr( glb_nb_nodes );
  array::ArrayView<uid_t,1> glb_id = array::make_view<uid_t,1>(glb_id_arr);

  parallel::mpi::comm().gatherv(loc_id.data(), loc_id.size(), glb_id.data(), recvcounts.data(), recvdispls.data(), root);
  // 2) Sort all global indices, and renumber from 1 to glb_nb_edges
  std::vector<Node> node_sort; node_sort.reserve(glb_nb_nodes);
  for( size_t jnode=0; jnode<glb_id.shape(0); ++jnode )
  {
    node_sort.push_back( Node(glb_id(jnode),jnode) );
  }
  std::sort(node_sort.begin(), node_sort.end());

  // Assume edge gid start
  uid_t gid=0;
  for( size_t jnode=0; jnode<node_sort.size(); ++jnode )
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

  parallel::mpi::comm().scatterv(glb_id.data(), recvcounts.data(), recvdispls.data(), loc_id.data(), loc_id.size(), root);

  for( int jnode=0; jnode<nb_nodes; ++jnode )
  {
    glb_idx(jnode) = loc_id(jnode);
  }
}

//----------------------------------------------------------------------------------------------------------------------

Field& build_nodes_remote_idx( mesh::Nodes& nodes )
{
  Timer scope_timer(__FUNCTION__);
  size_t mypart = parallel::mpi::comm().rank();
  size_t nparts = parallel::mpi::comm().size();

  UniqueLonLat compute_uid(nodes);

  // This piece should be somewhere central ... could be NPROMA ?
  // ---------->
  std::vector< int > proc( nparts );
  for( size_t jpart=0; jpart<nparts; ++jpart )
    proc[jpart] = jpart;
  // <---------

  array::IndexView<int,   1> ridx   = array::make_indexview<int,1>( nodes.remote_index()  );
  array::ArrayView<int,   1> part   = array::make_view<int,1>( nodes.partition()   );
  size_t nb_nodes = nodes.size();

  int varsize=2;

  std::vector< std::vector<uid_t> > send_needed( parallel::mpi::comm().size() );
  std::vector< std::vector<uid_t> > recv_needed( parallel::mpi::comm().size() );
  int sendcnt=0;
  std::map<uid_t,int> lookup;
  for( size_t jnode=0; jnode<nb_nodes; ++jnode )
  {
    uid_t uid = compute_uid(jnode);
    if( size_t(part(jnode)) == mypart )
    {
      lookup[ uid ] = jnode;
      ridx(jnode) = jnode;
    }
    else
    {
      ASSERT( jnode < part.shape(0) );
      if( part(jnode) >= (int)proc.size() )
      {
        std::stringstream msg;
        msg << "Assertion [part("<<jnode<<") < proc.size()] failed\n"
            << "part("<<jnode<<") = " << part(jnode) << "\n"
            << "proc.size() = " << proc.size();
        eckit::AssertionFailed(msg.str(),Here());
      }
      ASSERT( part(jnode) < (int)proc.size() );
      ASSERT( (size_t)proc[part(jnode)] < send_needed.size() );
      send_needed[ proc[part(jnode)] ].push_back( uid  );
      send_needed[ proc[part(jnode)] ].push_back( jnode );
      sendcnt++;
    }
  }

  parallel::mpi::comm().allToAll(send_needed, recv_needed);

  std::vector< std::vector<int> > send_found( parallel::mpi::comm().size() );
  std::vector< std::vector<int> > recv_found( parallel::mpi::comm().size() );

  for( size_t jpart=0; jpart<nparts; ++jpart )
  {
    const std::vector<uid_t>& recv_node = recv_needed[ proc[jpart] ];
    const size_t nb_recv_nodes = recv_node.size()/varsize;
    // array::ArrayView<uid_t,2> recv_node( make_view( Array::wrap(shape, recv_needed[ proc[jpart] ].data()) ),
    //     array::make_shape(recv_needed[ proc[jpart] ].size()/varsize,varsize) );
    for( size_t jnode=0; jnode<nb_recv_nodes; ++jnode )
    {
      uid_t uid = recv_node[jnode*varsize+0];
      int inode = recv_node[jnode*varsize+1];
      if( lookup.count(uid) )
      {
        send_found[ proc[jpart] ].push_back( inode );
        send_found[ proc[jpart] ].push_back( lookup[uid] );
      }
      else
      {
        std::stringstream msg;
        msg << "[" << parallel::mpi::comm().rank() << "] " << "Node requested by rank ["<<jpart
            << "] with uid [" << uid << "] that should be owned is not found";
        throw eckit::SeriousBug(msg.str(),Here());
      }
    }
  }

  parallel::mpi::comm().allToAll(send_found, recv_found);

  for( size_t jpart=0; jpart<nparts; ++jpart )
  {
    const std::vector<int>& recv_node = recv_found[ proc[jpart] ];
    const size_t nb_recv_nodes = recv_node.size()/2;
    // array::ArrayView<int,2> recv_node( recv_found[ proc[jpart] ].data(),
    //     array::make_shape(recv_found[ proc[jpart] ].size()/2,2) );
    for( size_t jnode=0; jnode<nb_recv_nodes; ++jnode )
    {
      ridx( recv_node[jnode*2+0] ) = recv_node[jnode*2+1];
    }
  }
  return nodes.remote_index();
}

//----------------------------------------------------------------------------------------------------------------------

Field& build_nodes_partition( mesh::Nodes& nodes )
{
  Timer scope_timer(__FUNCTION__);
  return nodes.partition();
}

//----------------------------------------------------------------------------------------------------------------------

Field& build_edges_partition( Mesh& mesh )
{
  Timer scope_timer(__FUNCTION__);

  const mesh::Nodes& nodes = mesh.nodes();

  UniqueLonLat compute_uid(mesh);

  size_t mypart = parallel::mpi::comm().rank();
  size_t nparts = parallel::mpi::comm().size();

  mesh::HybridElements& edges = mesh.edges();
  array::ArrayView<int,1> edge_part = array::make_view<int,1>( edges.partition() );
  const mesh::HybridElements::Connectivity& edge_nodes = edges.node_connectivity();
  const mesh::HybridElements::Connectivity& edge_to_elem = edges.cell_connectivity();

  std::shared_ptr< array::ArrayView<int,1> > is_pole_edge;
  bool has_pole_edges = false;
  if( edges.has_field("is_pole_edge") )
  {
    has_pole_edges = true;
    is_pole_edge = std::shared_ptr< array::ArrayView<int,1> > (
      new array::ArrayView<int,1>( array::make_view<int,1>( edges.field("is_pole_edge") ) ) );
  }

  array::ArrayView<int,1> node_part = array::make_view<int,1>( nodes.partition() );
  array::ArrayView<double,2> xy = array::make_view<double,2>( nodes.xy() );
  array::ArrayView<int,   1> flags  = array::make_view<int,1>( nodes.field("flags") );
#ifdef DEBUGGING_PARFIELDS
  array::ArrayView<gidx_t,   1> gidx  = array::make_view<gidx_t,1>( nodes.global_index() );
#endif

  array::ArrayView<int,1>     elem_part = array::make_view<int,1>(mesh.cells().partition() );

  PeriodicTransform transform;

  size_t nb_edges = edges.size();

  std::vector<int> periodic(nb_edges);

  for( size_t jedge=0; jedge<nb_edges; ++jedge )
  {
    periodic[jedge] = 0;
    idx_t ip1 = edge_nodes(jedge,0);
    idx_t ip2 = edge_nodes(jedge,1);
    int pn1 = node_part( ip1 );
    int pn2 = node_part( ip2 );
    if( pn1 == pn2 )
      edge_part(jedge) = pn1;
    else
    {
      if( edge_to_elem(jedge,1) == edge_to_elem.missing_value() ) // This is a edge at partition boundary
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
          int pe1 = elem_part(edge_to_elem(jedge,0));
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
        int pe1 = elem_part(edge_to_elem(jedge,0));
        int pe2 = elem_part(edge_to_elem(jedge,1));
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
    std::vector< std::vector<uid_t> > send_unknown( parallel::mpi::comm().size() );
    std::vector< std::vector<uid_t> > recv_unknown( parallel::mpi::comm().size() );
    for( size_t jedge=0; jedge<nb_edges; ++jedge )
    {
      int ip1 = edge_nodes(jedge,0);
      int ip2 = edge_nodes(jedge,1);
      int pn1 = node_part( ip1 );
      int pn2 = node_part( ip2 );

      centroid[XX] = 0.5*(xy( ip1, XX ) + xy( ip2, XX ) );
      centroid[YY] = 0.5*(xy( ip1, YY ) + xy( ip2, YY ) );
      if( has_pole_edges && (*is_pole_edge)(jedge) )
      {
        centroid[YY] = centroid[YY] > 0 ? 90. : -90.;
      }

      transform(centroid,periodic[jedge]);
      uid_t uid = util::unique_lonlat(centroid);

      if( size_t(edge_part(jedge)) == mypart )
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
          x1 = xy(ip1,XX);
          y1 = xy(ip1,YY);
          x2 = xy(ip2,XX);
          y2 = xy(ip2,YY);
          xe = centroid[XX];
          ye = centroid[YY];
          DEBUG( uid << " --> " << EDGE(jedge) << "   x1,y1 - x2,y2 - xe,ye " << x1<<","<<y1
                 << " - " << x2<<","<<y2<< " - " << xe <<","<<ye<< "     part " << edge_part(jedge));
        }
#endif

    }

    parallel::mpi::comm().allToAll(send_unknown, recv_unknown);

    // So now we have identified all possible edges with wrong partition.
    // We still need to check if it is actually wrong. This can be achieved
    // just by looking if the edge in question is owned by the possibly wrongly
    // assigned edge partition. If the edge is not found on that partition,
    // then its other node must be the edge partition.

    std::vector< std::vector<int> > send_found( parallel::mpi::comm().size() );
    std::vector< std::vector<int> > recv_found( parallel::mpi::comm().size() );

    for( size_t jpart=0; jpart<nparts; ++jpart )
    {
      const std::vector<uid_t>& recv_edge = recv_unknown[jpart];
      const size_t nb_recv_edges = recv_edge.size()/varsize;
      // array::ArrayView<uid_t,2> recv_edge( recv_unknown[ jpart ].data(),
      //     array::make_shape(recv_unknown[ jpart ].size()/varsize,varsize) );
      for( size_t jedge=0; jedge<nb_recv_edges; ++jedge )
      {
        uid_t uid      = recv_edge[jedge*varsize+0];
        int    recv_idx = recv_edge[jedge*varsize+1];
        if( lookup.count(uid) )
        {
          send_found[ jpart ].push_back( recv_idx );
        }
      }
    }

    parallel::mpi::comm().allToAll(send_found, recv_found);

    for( size_t jpart=0; jpart<nparts; ++jpart )
    {
      for( size_t jedge=0; jedge<recv_found[jpart].size(); ++jedge )
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

  // Sanity check
  for( size_t jedge=0; jedge<nb_edges; ++jedge )
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

  return edges.partition();
}


Field& build_edges_remote_idx( Mesh& mesh  )
{
  Timer scope_timer(__FUNCTION__);

  const mesh::Nodes& nodes = mesh.nodes();
  UniqueLonLat compute_uid(mesh);

  size_t mypart = parallel::mpi::comm().rank();
  size_t nparts = parallel::mpi::comm().size();

  mesh::HybridElements& edges = mesh.edges();

  array::IndexView<int,1>       edge_ridx = array::make_indexview<int,1> ( edges.remote_index() );

  const array::ArrayView<int,1> edge_part = array::make_view<int,1>( edges.partition() );
  const mesh::HybridElements::Connectivity& edge_nodes = edges.node_connectivity();

  array::ArrayView<double,2> xy = array::make_view<double,2>( nodes.xy() );
  array::ArrayView<int,   1> flags  = array::make_view<int,1>( nodes.field("flags")       );
#ifdef DEBUGGING_PARFIELDS
  array::ArrayView<gidx_t,   1> gidx   = array::make_view<gidx_t,1>( nodes.global_index()       );
  array::ArrayView<int,   1> node_part = array::make_view<int,1>( nodes.partition()       );
#endif

  std::shared_ptr< array::ArrayView<int,1> > is_pole_edge;
  bool has_pole_edges = false;
  if( edges.has_field("is_pole_edge") )
  {
    has_pole_edges = true;
    is_pole_edge = std::shared_ptr< array::ArrayView<int,1> > (
      new array::ArrayView<int,1>( array::make_view<int,1>( edges.field("is_pole_edge") ) ) );
  }

  const int nb_edges = edges.size();

  double centroid[2];
  std::vector< std::vector<uid_t> > send_needed( parallel::mpi::comm().size() );
  std::vector< std::vector<uid_t> > recv_needed( parallel::mpi::comm().size() );
  int sendcnt=0;
  std::map<uid_t,int> lookup;

  PeriodicTransform transform;


  for( int jedge=0; jedge<nb_edges; ++jedge )
  {
    int ip1 = edge_nodes(jedge,0);
    int ip2 = edge_nodes(jedge,1);
    centroid[XX] = 0.5*(xy( ip1, XX ) + xy( ip2, XX ) );
    centroid[YY] = 0.5*(xy( ip1, YY ) + xy( ip2, YY ) );
    if( has_pole_edges && (*is_pole_edge)(jedge) )
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

    uid_t uid = util::unique_lonlat(centroid);
    if( size_t(edge_part(jedge)) == mypart && !needed ) // All interior edges fall here
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
  parallel::mpi::comm().allToAll(send_needed, recv_needed);

  std::vector< std::vector<int> > send_found( parallel::mpi::comm().size() );
  std::vector< std::vector<int> > recv_found( parallel::mpi::comm().size() );

  std::map<uid_t,int>::iterator found;
  for( size_t jpart=0; jpart<nparts; ++jpart )
  {
    const std::vector<uid_t>& recv_edge = recv_needed[jpart];
    const size_t nb_recv_edges = recv_edge.size()/varsize;
    // array::ArrayView<uid_t,2> recv_edge( recv_needed[ jpart ].data(),
    //     array::make_shape(recv_needed[ jpart ].size()/varsize,varsize) );
    for( size_t jedge=0; jedge<nb_recv_edges; ++jedge )
    {
      uid_t recv_uid = recv_edge[jedge*varsize+0];
      int recv_idx = recv_edge[jedge*varsize+1];
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
        Log::warning() << msg.str() << " @ " << Here() << std::endl;
      }
    }
  }
  parallel::mpi::comm().allToAll(send_found, recv_found);
  for( size_t jpart=0; jpart<nparts; ++jpart )
  {
    const std::vector<int>& recv_edge = recv_found[jpart];
    const size_t nb_recv_edges = recv_edge.size()/2;
    // array::ArrayView<int,2> recv_edge( recv_found[ jpart ].data(),
    //     array::make_shape(recv_found[ jpart ].size()/2,2) );
    for( size_t jedge=0; jedge<nb_recv_edges; ++jedge )
    {
      edge_ridx( recv_edge[jedge*2+0] ) = recv_edge[jedge*2+1];
    }
  }
  return edges.remote_index();
}

Field& build_edges_global_idx( Mesh& mesh )
{
  Timer scope_timer(__FUNCTION__);

  UniqueLonLat compute_uid(mesh);

  int nparts = parallel::mpi::comm().size();
  size_t root = 0;

  mesh::HybridElements& edges = mesh.edges();

  array::make_storageview<gidx_t>(edges.global_index()).assign(-1);
  array::ArrayView<gidx_t,1> edge_gidx = array::make_view<gidx_t,1>( edges.global_index() );

  const mesh::HybridElements::Connectivity& edge_nodes = edges.node_connectivity();
  array::ArrayView<double,2> xy = array::make_view<double,2>( mesh.nodes().xy() );
  std::shared_ptr< array::ArrayView<int,1> > is_pole_edge;
  bool has_pole_edges = false;
  if( edges.has_field("is_pole_edge") )
  {
    has_pole_edges = true;
    is_pole_edge = std::shared_ptr< array::ArrayView<int,1> > (
      new array::ArrayView<int,1>( array::make_view<int,1>( edges.field("is_pole_edge") ) ) );
  }

  /*
   * Sorting following edge_gidx will define global order of
   * gathered fields. Special care needs to be taken for
   * pole edges, as their centroid might coincide with
   * other edges
   */
  double centroid[2];
  int nb_edges = edges.size();
  for( int jedge=0; jedge<nb_edges; ++jedge )
  {
    if( edge_gidx(jedge) <= 0 )
    {
      centroid[XX] = 0.5*(xy( edge_nodes(jedge,0), XX ) + xy( edge_nodes(jedge,1), XX ) );
      centroid[YY] = 0.5*(xy( edge_nodes(jedge,0), YY ) + xy( edge_nodes(jedge,1), YY ) );
      if( has_pole_edges && (*is_pole_edge)(jedge) )
      {
        centroid[YY] = centroid[YY] > 0 ? 90. : -90.;
      }
      edge_gidx(jedge) = util::unique_lonlat(centroid);
    }
  }

  /*
   * REMOTE INDEX BASE = 1
   */
  // unused //  const int ridx_base = 1;

  // 1) Gather all global indices, together with location
  array::ArrayT<uid_t> loc_edge_id_arr(nb_edges);
  array::ArrayView<uid_t,1> loc_edge_id = array::make_view<uid_t,1>(loc_edge_id_arr);

  for( int jedge=0; jedge<nb_edges; ++jedge )
  {
    loc_edge_id(jedge) = edge_gidx(jedge);
  }

  std::vector<int> recvcounts(parallel::mpi::comm().size());
  std::vector<int> recvdispls(parallel::mpi::comm().size());

  parallel::mpi::comm().gather(nb_edges, recvcounts, root);

  recvdispls[0]=0;
  for (int jpart=1; jpart<nparts; ++jpart) // start at 1
  {
    recvdispls[jpart]=recvcounts[jpart-1]+recvdispls[jpart-1];
  }
  int glb_nb_edges = std::accumulate(recvcounts.begin(),recvcounts.end(),0);

  array::ArrayT<uid_t> glb_edge_id_arr(glb_nb_edges);
  array::ArrayView<uid_t,1> glb_edge_id = array::make_view<uid_t,1>(glb_edge_id_arr);

  parallel::mpi::comm().gatherv(loc_edge_id.data(), loc_edge_id.size(), glb_edge_id.data(), recvcounts.data(), recvdispls.data(), root);
  // 2) Sort all global indices, and renumber from 1 to glb_nb_edges
  std::vector<Node> edge_sort; edge_sort.reserve(glb_nb_edges);
  for( size_t jedge=0; jedge<glb_edge_id.shape(0); ++jedge )
  {
    edge_sort.push_back( Node(glb_edge_id(jedge),jedge) );
  }
  std::sort(edge_sort.begin(), edge_sort.end());

  // Assume edge gid start
  uid_t gid=200000;
  for( size_t jedge=0; jedge<edge_sort.size(); ++jedge )
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

  parallel::mpi::comm().scatterv(glb_edge_id.data(), recvcounts.data(), recvdispls.data(), loc_edge_id.data(), loc_edge_id.size(), root);

  for( int jedge=0; jedge<nb_edges; ++jedge )
  {
    edge_gidx(jedge) = loc_edge_id(jedge);
  }

  return edges.global_index();
}

//----------------------------------------------------------------------------------------------------------------------
// C wrapper interfaces to C++ routines

void atlas__build_parallel_fields ( Mesh::Implementation* mesh) {
  ATLAS_ERROR_HANDLING( Mesh m(mesh); build_parallel_fields(m); );
}
void atlas__build_nodes_parallel_fields (mesh::Nodes* nodes) {
  ATLAS_ERROR_HANDLING( build_nodes_parallel_fields(*nodes) );
}

void atlas__build_edges_parallel_fields ( Mesh::Implementation* mesh ) {
  ATLAS_ERROR_HANDLING( Mesh m(mesh); build_edges_parallel_fields(m); );
}

void atlas__renumber_nodes_glb_idx (mesh::Nodes* nodes)
{
  ATLAS_ERROR_HANDLING( renumber_nodes_glb_idx(*nodes) );
}

//----------------------------------------------------------------------------------------------------------------------

} // namespace actions
} // namespace mesh
} // namespace atlas

