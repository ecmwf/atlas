/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */


/// @warning Still doesn't know about periodic BC to enlarge Halo

#include <iostream>
#include <stdexcept>
#include <cmath>
#include <limits>

#include "eckit/memory/ScopedPtr.h"

#include "atlas/mpi/mpi.h"
#include "atlas/atlas_config.h"
#include "atlas/runtime/ErrorHandling.h"
#include "atlas/Mesh.h"
#include "atlas/FunctionSpace.h"
#include "atlas/Field.h"
#include "atlas/Nodes.h"
#include "atlas/actions/BuildHalo.h"
#include "atlas/Parameters.h"
#include "atlas/util/ArrayView.h"
#include "atlas/util/IndexView.h"
#include "atlas/Array.h"
#include "atlas/util/AccumulateFaces.h"
#include "atlas/util/Bitflags.h"
#include "atlas/util/PeriodicTransform.h"
#include "atlas/util/Unique.h"
#include "atlas/util/LonLatMicroDeg.h"
#include "atlas/util/Functions.h"

using atlas::util::Face;
using atlas::util::accumulate_faces;
using atlas::util::Topology;
using atlas::util::PeriodicTransform;
using atlas::util::UniqueLonLat;
using atlas::util::LonLatMicroDeg;
using atlas::util::microdeg;

namespace atlas {
namespace actions {


typedef gidx_t uid_t;

// ------------------------------------------------------------------
class BuildHaloHelper;

void increase_halo( Mesh& mesh );
void increase_halo_interior( BuildHaloHelper& );




class EastWest: public PeriodicTransform
{
public:
	EastWest()
	{
		x_translation_ = -360.;
	}
};


class WestEast: public PeriodicTransform
{
public:
	WestEast()
	{
		x_translation_ = 360.;
	}
};


void increase_halo( Mesh& mesh )
{
  //DEBUG( "\n\n" << "Increase halo!! \n\n");
  Nodes& nodes         = mesh.nodes();
  ArrayView<double,2> lonlat   ( nodes.lonlat() );
  ArrayView<gidx_t,1> glb_idx  ( nodes.global_index() );
  ArrayView<int   ,1> part     ( nodes.partition() );
  IndexView<int   ,1> ridx     ( nodes.remote_index() );
  ArrayView<int   ,1> flags    ( nodes.field( "flags"          ) );

  size_t nb_nodes = nodes.size();

  PeriodicTransform transform;

  std::vector< std::vector< ElementRef > > node_to_elem(nb_nodes);

  std::vector< IndexView<int,2> > elem_nodes( mesh.nb_function_spaces() );
  std::vector< ArrayView<int,1> > elem_part ( mesh.nb_function_spaces() );
  std::vector< ArrayView<gidx_t,1> > elem_glb_idx ( mesh.nb_function_spaces() );

  for(size_t func_space_idx = 0; func_space_idx < mesh.nb_function_spaces(); ++func_space_idx)
  {
    FunctionSpace& elements = mesh.function_space(func_space_idx);
    if( elements.metadata().get<long>("type") == Entity::ELEMS )
    {
      elem_nodes  [func_space_idx] = IndexView<int,2>( elements.field("nodes") );
      elem_part   [func_space_idx] = ArrayView<int,1>( elements.field("partition") );
      elem_glb_idx[func_space_idx] = ArrayView<gidx_t,1>( elements.field("glb_idx") );
      int nb_elems = elem_nodes[func_space_idx].shape(0);
      int nb_nodes_per_elem = elem_nodes[func_space_idx].shape(1);
      for (int elem=0; elem<nb_elems; ++elem)
      {
        for (int n=0; n<nb_nodes_per_elem; ++n)
        {
          int node = elem_nodes[func_space_idx](elem,n);
          node_to_elem[node].push_back( ElementRef(elements.index(),elem) );
        }
      }
    }
  }

  UniqueLonLat compute_uid(nodes);


  /*
  1) Find nodes at boundary of partition
  2) Communicate glb_index of these boundary nodes to other partitions
  3) Find received glb_index in glb_node_to_local_node list
  4) Find elements in node_to_elem list that belong to me
  5) Make list of all nodes that complete the elements
  6) Communicate elements and nodes back
  7) Adapt mesh
  */


  /*
  1) Find boundary of partition:
  - find unique edges
  - if edge is bdry_edge, then the nodes are bdry nodes
  */
  std::set<int> bdry_nodes_set;

  FunctionSpace& quads       = mesh.function_space( "quads" );
  FunctionSpace& triags      = mesh.function_space( "triags" );

  std::vector< std::vector<int> > node_to_face(nb_nodes);
  std::vector< int > face_nodes_data; face_nodes_data.reserve(4*nb_nodes);
  std::vector< Face > face_to_elem;
  face_to_elem.reserve(4*nb_nodes);
  int nb_faces = 0;
  int nb_inner_faces = 0;

  accumulate_faces(quads, node_to_face,face_nodes_data,face_to_elem,nb_faces,nb_inner_faces);
  accumulate_faces(triags,node_to_face,face_nodes_data,face_to_elem,nb_faces,nb_inner_faces);

  size_t extents[] = {size_t(nb_faces), 2};
  ArrayView<int,2> face_nodes(face_nodes_data.data(),extents);


  for( int jface=0; jface<nb_faces; ++jface )
  {
    if( face_to_elem[jface].is_bdry() )
    {
      for( int jnode=0; jnode<2; ++jnode) // 2 nodes per face
      {
        if( face_nodes(jface,jnode) >= 0 )
        {
          bdry_nodes_set.insert(face_nodes(jface,jnode));
        }
      }
    }
  }

  std::vector<int> bdry_nodes( bdry_nodes_set.begin(), bdry_nodes_set.end());

  /*
  2) Communicate glb_index of these boundary nodes to other partitions
  3) Find received glb_index in glb_node_to_local_node list
  */

  std::map<uid_t,int> node_uid_to_loc;
  for( size_t jnode=0; jnode<nb_nodes; ++jnode )
  {
    LonLatMicroDeg ll(lonlat[jnode]);
    if( node_uid_to_loc.count(ll.unique()) > 0 )
    {
      int other = node_uid_to_loc[ll.unique()];
      std::stringstream msg;
      msg << "Node uid: " << ll.unique() << "   " << glb_idx(jnode)
          << " (" << lonlat(jnode,LON) <<","<< lonlat(jnode,LAT)<<")  has already been added as node "
          << glb_idx(other) << " (" << lonlat(other,LON) <<","<< lonlat(other,LAT)<<")";
      throw eckit::SeriousBug(msg.str(),Here());
    }
    node_uid_to_loc[ll.unique()] = jnode;
  }

  size_t nb_bdry_nodes = bdry_nodes.size();
  eckit::ScopedPtr<Array> arr_bdry_nodes_id( Array::create<uid_t>(nb_bdry_nodes,4) );
  ArrayView<uid_t,2> bdry_nodes_id(*arr_bdry_nodes_id);
  ASSERT( bdry_nodes_id.shape(0) == nb_bdry_nodes );
  ASSERT( bdry_nodes_id.shape(1) == 4);

  for(size_t jnode = 0; jnode < nb_bdry_nodes; ++jnode)
  {
    LonLatMicroDeg ll( lonlat[bdry_nodes[jnode]] );
    bdry_nodes_id(jnode,0) = ll.lon();
    bdry_nodes_id(jnode,1) = ll.lat();
    bdry_nodes_id(jnode,2) = glb_idx( bdry_nodes[jnode] );
    bdry_nodes_id(jnode,3) = flags( bdry_nodes[jnode] );
  }

  std::vector<int> recvcounts( eckit::mpi::size() );
  std::vector<int> recvdispls( eckit::mpi::size() );
  int sendcnt = bdry_nodes_id.size();
  ASSERT( size_t(sendcnt) == nb_bdry_nodes*4 );
  ECKIT_MPI_CHECK_RESULT( MPI_Allgather( &sendcnt,          1, MPI_INT,
                                   recvcounts.data(), 1, MPI_INT, eckit::mpi::comm() ) );

  recvdispls[0] = 0;
  int recvcnt = recvcounts[0];
  for( size_t jpart=1; jpart<eckit::mpi::size(); ++jpart )
  {
    recvdispls[jpart] = recvdispls[jpart-1] + recvcounts[jpart-1];
    recvcnt += recvcounts[jpart];
  }
  std::vector<uid_t> recvbuf(recvcnt);

  ECKIT_MPI_CHECK_RESULT( MPI_Allgatherv( bdry_nodes_id.data(), sendcnt, eckit::mpi::datatype<uid_t>(),
                    recvbuf.data(), recvcounts.data(), recvdispls.data(),
                    eckit::mpi::datatype<uid_t>(), eckit::mpi::comm()) );

  // sfn stands for "send_found_nodes"
  std::vector< std::vector<int>    > sfn_part( eckit::mpi::size() );
  std::vector< std::vector<int>    > sfn_ridx( eckit::mpi::size() );
  std::vector< std::vector<uid_t> > sfn_glb_idx( eckit::mpi::size() );
  std::vector< std::vector<int>    > sfn_flags ( eckit::mpi::size() );
  std::vector< std::vector<double> > sfn_lonlat ( eckit::mpi::size() );
  // sfn stands for "send_found_elems"
  std::vector< std::vector< std::vector<uid_t> > >
      sfe_glb_idx ( mesh.nb_function_spaces(), std::vector< std::vector<uid_t> >( eckit::mpi::size() ) );
  std::vector< std::vector< std::vector<uid_t> > >
      sfe_nodes_id( mesh.nb_function_spaces(), std::vector< std::vector<uid_t> >( eckit::mpi::size() ) );
  std::vector< std::vector< std::vector<int> > >
      sfe_part    ( mesh.nb_function_spaces(), std::vector< std::vector<int> >( eckit::mpi::size() ) );

  // 4) Find elements in node_to_elem list that belong to me
  // 5) Make list of all nodes that complete the elements

  for (size_t jpart=0; jpart<eckit::mpi::size(); ++jpart)
  {
    ArrayView<uid_t,2> recv_bdry_nodes_id( recvbuf.data()+recvdispls[jpart],
                                         make_shape( recvcounts[jpart]/4, 4 ).data() );
    size_t recv_nb_bdry_nodes = recv_bdry_nodes_id.shape(0);

    // Find elements that have these nodes
    // In order to do this, check the node_to_elem list
    // Warning: only add elements that are owned!


    std::vector< std::set< std::pair<int,int> > > found_bdry_elements_set( mesh.nb_function_spaces() );
    for(size_t jrecv = 0; jrecv < recv_nb_bdry_nodes; ++jrecv)
    {
      ASSERT( recv_bdry_nodes_id.shape(1) == 4 );
      int recv_x       = recv_bdry_nodes_id(jrecv,0);
      int recv_y       = recv_bdry_nodes_id(jrecv,1);
      uid_t recv_glb_idx = recv_bdry_nodes_id(jrecv,2);
      int recv_flags   = recv_bdry_nodes_id(jrecv,3);

      LonLatMicroDeg ll(recv_x,recv_y);

      int periodic=0;
      // If the received node is flagged as periodic, look for point following periodic transformation
      if( Topology::check(recv_flags,Topology::PERIODIC|Topology::EAST) ) periodic = -1; // If the node is at BC_EAST (so slave), the master is in negative direction (-360 deg)
      else if( Topology::check(recv_flags,Topology::PERIODIC|Topology::WEST) )
        periodic =  1; // If the node is at BC_WEST (so master), the slave is in positive direction (+360 deg)

      std::vector<uid_t> recv_uids; recv_uids.reserve(2);
      recv_uids.push_back( ll.unique() );
      if( periodic ) {
        transform( ll, periodic );
        recv_uids.push_back( ll.unique() );
      }

      for( size_t juid = 0; juid < recv_uids.size(); ++juid )
      {
        uid_t recv_uid = recv_uids[juid];
        int loc = -1;

        // search and get local node index for received node
        std::map<uid_t,int>::iterator found = node_uid_to_loc.find(recv_uid);
        if( found != node_uid_to_loc.end() )
        {
          loc = found->second;
          if( eckit::mpi::rank() == jpart && glb_idx(loc) == recv_glb_idx ) loc = -1;
        }
        //if( periodic && loc != -1 ) DEBUG(" found it at " << glb_idx(loc));

        if( loc != -1 )
        {
          for( size_t jelem=0; jelem<node_to_elem[loc].size(); ++jelem )
          {
            int f = node_to_elem[loc][jelem].f;
            int e = node_to_elem[loc][jelem].e;
            if( size_t(elem_part[f](e)) == eckit::mpi::rank() )
            {
              //DEBUG( "node " << recv_glb_idx << "\t  --> " << elem_glb_idx[f][e] );
              found_bdry_elements_set[f].insert( std::make_pair( e, (juid==0 ? 0 : -periodic) ) );
              //if( periodic && loc != -1 ) DEBUG(" going to send element " << elem_glb_idx[f][e]);
            }
          }
        }
      }
    }
    // found_bdry_elements_set now contains elements for the nodes

    std::vector< std::vector<int> > found_bdry_elements( mesh.nb_function_spaces() );
    std::vector< std::vector<int> > found_bdry_elements_periodic( mesh.nb_function_spaces() );

    std::vector<int> nb_found_bdry_elems( mesh.nb_function_spaces(), 0 );

    for(size_t f = 0; f < mesh.nb_function_spaces(); ++f)
    {
      nb_found_bdry_elems[f] = found_bdry_elements_set[f].size();

      found_bdry_elements[f].resize(nb_found_bdry_elems[f]);
      found_bdry_elements_periodic[f].resize(nb_found_bdry_elems[f]);
      int jelem=0;
      std::set< std::pair<int,int> >::iterator it=found_bdry_elements_set[f].begin();
      for( ; it!=found_bdry_elements_set[f].end(); ++it, ++jelem )
      {
        found_bdry_elements[f][jelem] = it->first;
        found_bdry_elements_periodic[f][jelem] = it->second;
      }
      nb_found_bdry_elems[f] = found_bdry_elements[f].size();
    }

    // Collect all nodes needed to complete the element, and mark if they become periodic on requesting task
    std::set< std::pair<LonLatMicroDeg,int> > found_bdry_nodes_id_set;
    {
      for( size_t f=0; f<mesh.nb_function_spaces(); ++f )
      {
        for( int jelem=0; jelem<nb_found_bdry_elems[f]; ++jelem )
        {
          int e = found_bdry_elements[f][jelem];
          int nb_elem_nodes = elem_nodes[f].shape(1);
          for( int n=0; n<nb_elem_nodes; ++n )
          {
            int x = microdeg( lonlat( elem_nodes[f](e,n), LON) );
            int y = microdeg( lonlat( elem_nodes[f](e,n), LAT) );
            int periodic = found_bdry_elements_periodic[f][jelem];
            found_bdry_nodes_id_set.insert( std::make_pair(LonLatMicroDeg(x,y),periodic) );
          }
        }
      }

      // Remove nodes we already have received, as we won't need to send them
      for( size_t jrecv=0; jrecv<recv_nb_bdry_nodes; ++jrecv)
      {
        int x = recv_bdry_nodes_id(jrecv,0);
        int y = recv_bdry_nodes_id(jrecv,1);
        int periodic=0;
        if( Topology::check(recv_bdry_nodes_id(jrecv,3),Topology::PERIODIC|Topology::EAST) )
          periodic = -1; // If the node is ghost (so slave), the master is in negative direction (-360 deg)
        else if( Topology::check(recv_bdry_nodes_id(jrecv,3),Topology::PERIODIC|Topology::WEST) )
          periodic = 1; // If the node is not ghost (so master), the slave is in positive direction (+360 deg)

        LonLatMicroDeg ll(x,y);
        found_bdry_nodes_id_set.erase( std::make_pair(ll,periodic) ) ;
        // DO I HAVE TO ALSO CHECK FOR PERIODICITY HERE?
      }
    }
    int nb_found_bdry_nodes = found_bdry_nodes_id_set.size();
    sfn_glb_idx[jpart].resize(nb_found_bdry_nodes);
    sfn_part[jpart].resize(nb_found_bdry_nodes);
    sfn_ridx[jpart].resize(nb_found_bdry_nodes);
    sfn_flags[jpart].resize(nb_found_bdry_nodes,Topology::NONE);
    sfn_lonlat[jpart].resize(2*nb_found_bdry_nodes);
    //DEBUG_VAR( nb_found_bdry_nodes );
    // Fill buffers to send
    {
      int jnode=0;
      std::set<std::pair<LonLatMicroDeg,int> >::iterator it;
      for( it=found_bdry_nodes_id_set.begin(); it!=found_bdry_nodes_id_set.end(); ++it, ++jnode )
      {
        LonLatMicroDeg ll = it->first;
        int periodic = it->second;
        //eckit::Log::warning() << "\n" << "Looking for node with lonlat " << ll.x*1.e-6*180./M_PI << "," << ll.y*1.e-6*180./M_PI << ".  periodic = " << periodic << std::endl;

        //DEBUG_VAR( periodic );
        uid_t uid = ll.unique();

        std::map<uid_t,int>::iterator found = node_uid_to_loc.find( uid );
        if( found != node_uid_to_loc.end() ) // Point exists inside domain
        {
          int loc = found->second;
          sfn_glb_idx[jpart][jnode]      = glb_idx(loc);
          sfn_part   [jpart][jnode]      = part   (loc);
          sfn_ridx   [jpart][jnode]      = ridx   (loc);
          sfn_lonlat [jpart][jnode*2+LON] = lonlat (loc,LON);
          sfn_lonlat [jpart][jnode*2+LAT] = lonlat (loc,LAT);
          //DEBUG_VAR(glb_idx(loc));
          if( periodic )
          {
            if( periodic > 0 )
              Topology::set(sfn_flags[jpart][jnode],Topology::EAST);
            else
              Topology::set(sfn_flags[jpart][jnode],Topology::WEST);
            if( Topology::check(flags(loc),Topology::BC ) )
            {
              Topology::set(sfn_flags[jpart][jnode],Topology::BC );
              if( !Topology::check(flags(loc),Topology::EAST) )
                Topology::set(sfn_flags[jpart][jnode],Topology::PERIODIC);
            }
            else
            {
              Topology::set(sfn_flags[jpart][jnode],Topology::PERIODIC);
            }
            transform(&sfn_lonlat[jpart][2*jnode],(double) periodic);
            // The glb_idx is based on the destination location
            sfn_glb_idx[jpart][jnode] = util::unique_lonlat(sfn_lonlat[jpart][jnode*2+LON],sfn_lonlat[jpart][jnode*2+LAT]);
          }
          else
          {
						// When a node in a partition corner is found through an interior element, but touches the east boundary, it should
						// also be flagged as periodic and bc_east
						Topology::set(sfn_flags[jpart][jnode],flags(loc));
          }
        }
        else
        {
          eckit::Log::warning() << "Node needed by ["<<jpart<<"] with lonlat " << ll.lon()*1.e-6 << ","
                                << ll.lat()*1.e-6 << " was not found in ["<<eckit::mpi::rank()<<"]." << std::endl;
          ASSERT(false);
        }
//        if( periodic != 0 )
//        {
//           eckit::Log::warning() << "Node needed by ["<<jpart<<"] with lonlat " << ll.x*1.e-6 << "," << ll.y*1.e-6 << " was not found in ["<<eckit::mpi::rank()<<"]." << std::endl;
//        }
      }
    }

    for( size_t f=0; f<mesh.nb_function_spaces(); ++f )
    {
      FunctionSpace& elements = mesh.function_space(f);

      if( elements.metadata().get<long>("type") == Entity::ELEMS )
      {
        int nb_elem_nodes(elem_nodes[f].shape(1));
        sfe_glb_idx [f][jpart].resize( nb_found_bdry_elems[f] );
        sfe_part    [f][jpart].resize( nb_found_bdry_elems[f] );
        sfe_nodes_id[f][jpart].resize( nb_found_bdry_elems[f]*nb_elem_nodes );

        ArrayView<uid_t,2> sfe_nodes_id_view( sfe_nodes_id[f][jpart].data(),
                                            make_shape(nb_found_bdry_elems[f],nb_elem_nodes).data() );

        for( int jelem=0; jelem<nb_found_bdry_elems[f]; ++jelem )
        {
          int e = found_bdry_elements[f][jelem];
          int periodic = found_bdry_elements_periodic[f][jelem];
          sfe_part[f][jpart][jelem]    = elem_part[f][e];

          double centroid[2];
          centroid[LON] = 0.;
          centroid[LAT] = 0.;
          for( int n=0; n<nb_elem_nodes; ++n)
          {
            double crd[2];
            crd[LON] = lonlat(elem_nodes[f](e,n),LON);
            crd[LAT] = lonlat(elem_nodes[f](e,n),LAT);
            transform(crd, (double)periodic);
            centroid[LON] += crd[LON];
            centroid[LAT] += crd[LAT];
            sfe_nodes_id_view(jelem,n) = util::unique_lonlat(crd);
          }
          centroid[LON] /= static_cast<double>(nb_elem_nodes);
          centroid[LAT] /= static_cast<double>(nb_elem_nodes);
          sfe_glb_idx[f][jpart][jelem] = util::unique_lonlat(centroid) ;
        }
      }
    }
  }

  // 6) Now communicate all found fields back

  //    rfn stands for "recv_found_nodes"
  std::vector< std::vector<uid_t>  > rfn_glb_idx(eckit::mpi::size());
  std::vector< std::vector<int>    > rfn_part(eckit::mpi::size());
  std::vector< std::vector<int>    > rfn_ridx( eckit::mpi::size() );
  std::vector< std::vector<int>    > rfn_flags( eckit::mpi::size() );
  std::vector< std::vector<double> > rfn_lonlat(eckit::mpi::size());
  //    rfe stands for "recv_found_elems"
  std::vector< std::vector< std::vector<uid_t> > >
      rfe_glb_idx ( mesh.nb_function_spaces(), std::vector< std::vector<uid_t> >( eckit::mpi::size() ) );
  std::vector< std::vector< std::vector<int> > >
      rfe_part    ( mesh.nb_function_spaces(), std::vector< std::vector<int> >( eckit::mpi::size() ) );
  std::vector< std::vector< std::vector<uid_t> > >
      rfe_nodes_id( mesh.nb_function_spaces(), std::vector< std::vector<uid_t> >( eckit::mpi::size() ) );

  eckit::mpi::all_to_all(sfn_glb_idx,  rfn_glb_idx);
  eckit::mpi::all_to_all(sfn_part,     rfn_part);
  eckit::mpi::all_to_all(sfn_ridx,     rfn_ridx);
  eckit::mpi::all_to_all(sfn_flags,    rfn_flags);
  eckit::mpi::all_to_all(sfn_lonlat,   rfn_lonlat);
  for( size_t f=0; f<mesh.nb_function_spaces(); ++f )
  {
    eckit::mpi::all_to_all(sfe_glb_idx [f], rfe_glb_idx [f] );
    eckit::mpi::all_to_all(sfe_nodes_id[f], rfe_nodes_id[f] );
    eckit::mpi::all_to_all(sfe_part    [f], rfe_part    [f] );
  }


  // We now have everything we need in rfe_ and rfn_ vectors
  // Now adapt the mesh

  // Nodes might be duplicated from different Tasks. We need to identify unique entries
  std::set<uid_t> node_uid;
  for(size_t jnode = 0; jnode < nb_nodes; ++jnode)
  {
    node_uid.insert( util::unique_lonlat(lonlat[jnode]) );
  }
  std::vector< std::vector<int> > rfn_idx(eckit::mpi::size());
  for( size_t jpart=0; jpart<eckit::mpi::size(); ++jpart )
  {
    rfn_idx[jpart].reserve(rfn_glb_idx[jpart].size());
  }

  int nb_new_nodes=0;
  for( size_t jpart=0; jpart<eckit::mpi::size(); ++jpart )
  {
    for( size_t n=0; n<rfn_glb_idx[jpart].size(); ++n )
    {
      double x = rfn_lonlat[jpart][n*2+LON];
      double y = rfn_lonlat[jpart][n*2+LAT];
      bool inserted = node_uid.insert( util::unique_lonlat( x, y ) ).second;
      if( inserted ) {
        rfn_idx[jpart].push_back(n);
      }
    }
    nb_new_nodes += rfn_idx[jpart].size();
  }

  //DEBUG_VAR(nb_new_nodes);

  nodes.resize( nb_nodes+nb_new_nodes );

  flags   = ArrayView<int,   1>( nodes.field("flags") );
  glb_idx = ArrayView<gidx_t,1>( nodes.global_index() );
  part    = ArrayView<int,   1>( nodes.partition() );
  ridx    = IndexView<int,   1>( nodes.remote_index() );
  lonlat  = ArrayView<double,2>( nodes.lonlat() );


  int new_node=0;
  for( size_t jpart=0; jpart<eckit::mpi::size(); ++jpart )
  {
    for( size_t n=0; n<rfn_idx[jpart].size(); ++n )
    {
      int loc_idx = nb_nodes+new_node;
      Topology::reset(flags(loc_idx),rfn_flags[jpart][rfn_idx[jpart][n]]);
      glb_idx(loc_idx)    = rfn_glb_idx [jpart][rfn_idx[jpart][n]];
      part   (loc_idx)    = rfn_part    [jpart][rfn_idx[jpart][n]];
      ridx   (loc_idx)    = rfn_ridx    [jpart][rfn_idx[jpart][n]];
      lonlat (loc_idx,LON) = rfn_lonlat  [jpart][rfn_idx[jpart][n]*2+LON];
      lonlat (loc_idx,LAT) = rfn_lonlat  [jpart][rfn_idx[jpart][n]*2+LAT];
      uid_t uid = util::unique_lonlat( lonlat[loc_idx] );

      // make sure new node was not already there
      std::map<uid_t,int>::iterator found = node_uid_to_loc.find(uid);
      if( found != node_uid_to_loc.end() )
      {
        int other = found->second;
        std::stringstream msg;
        msg << "New node with uid " << uid << ":\n"  << glb_idx(loc_idx)
            << "("<<lonlat(loc_idx,LON)<<","<<lonlat(loc_idx,LAT)<<")\n";
        msg << "Existing already loc "<< other << "  :  " << glb_idx(other)
            << "("<<lonlat(other,LON)<<","<<lonlat(other,LAT)<<")\n";
        throw eckit::SeriousBug(msg.str(),Here());
      }
      node_uid_to_loc[ uid ] = nb_nodes+new_node;
      ++new_node;
    }
  }

  for( size_t f=0; f<mesh.nb_function_spaces(); ++f )
  {
    FunctionSpace& elements = mesh.function_space(f);
    if( elements.metadata().get<long>("type") == Entity::ELEMS )
    {

      std::set<uid_t> elem_uid;
      size_t nb_elems = elements.shape(0);
      for(size_t jelem = 0; jelem < nb_elems; ++jelem)
      {
        elem_uid.insert( compute_uid(elem_nodes[f][jelem]) );
      }

      std::vector< std::vector<int> > received_new_elems(eckit::mpi::size());
      for( size_t jpart=0; jpart<eckit::mpi::size(); ++jpart )
      {
        received_new_elems[jpart].reserve(rfe_glb_idx[f][jpart].size());
      }

      size_t nb_new_elems=0;
      for( size_t jpart=0; jpart<eckit::mpi::size(); ++jpart )
      {
        for( size_t e=0; e<rfe_glb_idx[f][jpart].size(); ++e )
        {
          bool inserted = elem_uid.insert( rfe_glb_idx[f][jpart][e] ).second;
          if( inserted )
          {
            received_new_elems[jpart].push_back(e);
          }
        }
        nb_new_elems += received_new_elems[jpart].size();
      }

      //DEBUG_VAR( nb_new_elems );

      int nb_nodes_per_elem = elem_nodes[f].shape(1);
      elements.resize( make_shape( nb_elems+nb_new_elems, FunctionSpace::UNDEF_VARS ) );
      elem_glb_idx[f] = ArrayView<gidx_t,1>( elements.field("glb_idx") );
      elem_nodes[f]   = IndexView<int,2>( elements.field("nodes")   );
      elem_part[f]    = ArrayView<int,1>( elements.field("partition")   );
      int new_elem=0;
      for( size_t jpart=0; jpart<eckit::mpi::size(); ++jpart )
      {
        for( size_t e=0; e<received_new_elems[jpart].size(); ++e )
        {
          int jelem = received_new_elems[jpart][e];
          elem_glb_idx[f](nb_elems+new_elem)   = rfe_glb_idx[f][jpart][jelem];
          elem_part   [f](nb_elems+new_elem)   = rfe_part[f][jpart][jelem];
          for( int n=0; n<nb_nodes_per_elem; ++n )
            elem_nodes[f](nb_elems+new_elem,n) = node_uid_to_loc[ rfe_nodes_id[f][jpart][jelem*nb_nodes_per_elem+n] ];
          ++new_elem;
        }
      }
    }
  }
}


typedef std::vector< std::vector< ElementRef > > Node2Elem;


void build_lookup_node2elem( const Mesh& mesh, Node2Elem& node2elem )
{
  const Nodes& nodes  = mesh.nodes();

  node2elem.resize(nodes.size());
  for( size_t jnode=0; jnode<node2elem.size(); ++jnode )
    node2elem[jnode].clear();

  std::vector< IndexView<int,2> > elem_nodes( mesh.nb_function_spaces() );

  for( size_t func_space_idx=0; func_space_idx<mesh.nb_function_spaces(); ++func_space_idx)
  {
    FunctionSpace& elements = mesh.function_space(func_space_idx);
    if( elements.metadata().get<long>("type") == Entity::ELEMS )
    {
      elem_nodes  [func_space_idx] = IndexView<int,2>( elements.field("nodes") );
      int nb_elems = elem_nodes[func_space_idx].shape(0);
      int nb_nodes_per_elem = elem_nodes[func_space_idx].shape(1);
      for (int elem=0; elem<nb_elems; ++elem)
      {
        for (int n=0; n<nb_nodes_per_elem; ++n)
        {
          int node = elem_nodes[func_space_idx](elem,n);
          node2elem[node].push_back( ElementRef(elements.index(),elem) );
        }
      }
    }
  }
}

struct NoFilter {
  bool operator()(int) const { return true; }
};

void accumulate_partition_bdry_nodes( Mesh& mesh, std::vector<int>& bdry_nodes )
{
  std::set<int> bdry_nodes_set;

  Nodes& nodes       = mesh.nodes();
  FunctionSpace& quads       = mesh.function_space( "quads" );
  FunctionSpace& triags      = mesh.function_space( "triags" );

  int nb_nodes = nodes.size();
  std::vector< std::vector<int> > node_to_face(nb_nodes);
  std::vector< int > face_nodes_data; face_nodes_data.reserve(4*nb_nodes);
  std::vector< Face > face_to_elem;
  face_to_elem.reserve(4*nb_nodes);
  int nb_faces = 0;
  int nb_inner_faces = 0;

  accumulate_faces(quads, node_to_face,face_nodes_data,face_to_elem,nb_faces,nb_inner_faces);
  accumulate_faces(triags,node_to_face,face_nodes_data,face_to_elem,nb_faces,nb_inner_faces);

  size_t extents[] = {size_t(nb_faces), 2};
  ArrayView<int,2> face_nodes(face_nodes_data.data(),extents);

  for( int jface=0; jface<nb_faces; ++jface )
  {
    if( face_to_elem[jface].is_bdry() )
    {
      for( int jnode=0; jnode<2; ++jnode) // 2 nodes per face
      {
        if( face_nodes(jface,jnode) >= 0 )
        {
          bdry_nodes_set.insert(face_nodes(jface,jnode));
        }
      }
    }
  }
  bdry_nodes = std::vector<int>( bdry_nodes_set.begin(), bdry_nodes_set.end());
}

template< typename Predicate >
std::vector<int> filter_nodes(std::vector<int> nodes, const Predicate& predicate )
{
  std::vector<int> filtered; filtered.reserve(nodes.size());
  for( int jnode=0; jnode<nodes.size(); ++jnode )
  {
    int inode = nodes[jnode];
    if( predicate(inode) )
      filtered.push_back(inode);
  }
  return filtered;
}

class Notification
{
public:
  void add_error(const std::string& note, const eckit::CodeLocation& loc )
  {
    notes.push_back( note + " @ " + std::string(loc) );
  }
  void add_error(const std::string& note )
  {
    notes.push_back( note );
  }

  bool error() const { return notes.size() > 0; }
  void reset() { notes.clear(); }

  std::string str() const
  {
    std::stringstream stream;
    for(size_t jnote = 0; jnote < notes.size(); ++jnote)
    {
      if ( jnote > 0 ) stream << "\n";
      stream << notes[jnote];
    }
    return stream.str();
  }

  operator std::string() const { return str(); }

private:
  friend std::ostream& operator<<(std::ostream& s, const Notification& notes) { s << notes.str();  return s; }

private:
  std::vector<std::string> notes;
};


typedef std::map<uid_t,int> Uid2Node;
void build_lookup_uid2node( Mesh& mesh, Uid2Node& uid2node )
{
  Notification notes;
  Nodes& nodes         = mesh.nodes();
  ArrayView<double,2> lonlat   ( nodes.lonlat() );
  ArrayView<gidx_t,1> glb_idx  ( nodes.global_index() );
  int nb_nodes = nodes.size();

  UniqueLonLat compute_uid(nodes);

  uid2node.clear();
  for( int jnode=0; jnode<nb_nodes; ++jnode )
  {
    uid_t uid = compute_uid(jnode);

    if( uid2node.count(uid) > 0 )
    {
      int other = uid2node[uid];
      std::stringstream msg;
      msg << "Node uid: " << uid << "   " << glb_idx(jnode)
          << " (" << lonlat(jnode,LON) <<","<< lonlat(jnode,LAT)<<")  has already been added as node "
          << glb_idx(other) << " (" << lonlat(other,LON) <<","<< lonlat(other,LAT)<<")";
      notes.add_error(msg.str());
    }
    uid2node[uid] = jnode;
  }
  if( notes.error() )
    throw eckit::SeriousBug(notes.str(),Here());
}

void accumulate_elements( const Mesh& mesh,
                          const ArrayView<uid_t,1>& node_uid,
                          const Uid2Node& uid2node,
                          const Node2Elem& node2elem,
                          std::vector< std::vector<int> >& found_elements,
                          std::set< uid_t >& new_nodes_uid
                        )
{
  std::vector< IndexView<int,2> > elem_nodes( mesh.nb_function_spaces() );
  std::vector< ArrayView<int,1> > elem_part ( mesh.nb_function_spaces() );

  for(size_t func_space_idx = 0; func_space_idx < mesh.nb_function_spaces(); ++func_space_idx)
  {
    FunctionSpace& elements = mesh.function_space(func_space_idx);
    if( elements.metadata().get<long>("type") == Entity::ELEMS )
    {
      elem_nodes  [func_space_idx] = IndexView<int,2>( elements.field("nodes") );
      elem_part   [func_space_idx] = ArrayView<int,1>( elements.field("partition") );
    }
  }

  int nb_nodes = node_uid.size();

  std::vector< std::set< int > > found_elements_set( mesh.nb_function_spaces() );

  for( int jnode=0; jnode<nb_nodes; ++jnode )
  {
    uid_t uid = node_uid(jnode);

    int inode = -1;
    // search and get node index for uid
    Uid2Node::const_iterator found = uid2node.find(uid);
    if( found != uid2node.end() )
    {
      inode = found->second;
    }
    if( inode != -1 && size_t(inode) < node2elem.size() )
    {
      for(size_t jelem = 0; jelem < node2elem[inode].size(); ++jelem)
      {
        int f = node2elem[inode][jelem].f;
        int e = node2elem[inode][jelem].e;
        if( size_t(elem_part[f](e)) == eckit::mpi::rank() )
        {
          found_elements_set[f].insert( e );
        }
      }
    }
  }

  // found_bdry_elements_set now contains elements for the nodes
  found_elements.resize(mesh.nb_function_spaces());
  for(size_t f = 0; f < mesh.nb_function_spaces(); ++f)
  {
    found_elements[f] = std::vector<int>( found_elements_set[f].begin(), found_elements_set[f].end());
  }

  UniqueLonLat compute_uid(mesh.nodes());

  // Collect all nodes
  new_nodes_uid.clear();
  for(size_t f = 0; f < mesh.nb_function_spaces(); ++f)
  {
    for(size_t jelem = 0; jelem < found_elements[f].size(); ++jelem)
    {
      int e = found_elements[f][jelem];
      int nb_elem_nodes = elem_nodes[f].shape(1);
      for( int n=0; n<nb_elem_nodes; ++n )
      {
        new_nodes_uid.insert( compute_uid(elem_nodes[f](e,n)));
      }
    }
  }

  // Remove nodes we already have
  for( int jnode=0; jnode<nb_nodes; ++jnode)
  {
    new_nodes_uid.erase( node_uid(jnode) ) ;
  }
}

class BuildHaloHelper
{
public:
  struct Buffers
  {
    std::vector< std::vector<int>    > node_part;

    std::vector< std::vector<int>    > node_ridx;

    std::vector< std::vector<int>    > node_flags;

    std::vector< std::vector<uid_t>  > node_glb_idx;

    std::vector< std::vector<double> > node_lonlat;

    std::vector< std::vector< std::vector<uid_t> > > elem_glb_idx;

    std::vector< std::vector< std::vector<uid_t> > > elem_nodes_id;

    std::vector< std::vector< std::vector<int> > > elem_part;

    Buffers(Mesh& mesh)
    {
      node_part.resize(eckit::mpi::size());
      node_ridx.resize(eckit::mpi::size());
      node_flags.resize(eckit::mpi::size());
      node_glb_idx.resize(eckit::mpi::size());
      node_lonlat.resize(eckit::mpi::size());
      elem_glb_idx.resize( mesh.nb_function_spaces(), std::vector< std::vector<uid_t> >( eckit::mpi::size() ) );
      elem_nodes_id.resize( mesh.nb_function_spaces(), std::vector< std::vector<uid_t> >( eckit::mpi::size() ) );
      elem_part.resize( mesh.nb_function_spaces(), std::vector< std::vector<int> >( eckit::mpi::size() ) );
    }
  };
  static void all_to_all(Buffers& send, Buffers& recv)
  {
    eckit::mpi::all_to_all(send.node_glb_idx,  recv.node_glb_idx);
    eckit::mpi::all_to_all(send.node_part,     recv.node_part);
    eckit::mpi::all_to_all(send.node_ridx,     recv.node_ridx);
    eckit::mpi::all_to_all(send.node_flags,    recv.node_flags);
    eckit::mpi::all_to_all(send.node_lonlat,   recv.node_lonlat);
    for(size_t f = 0; f < send.elem_glb_idx.size(); ++f)
    {
      eckit::mpi::all_to_all(send.elem_glb_idx [f], recv.elem_glb_idx [f] );
      eckit::mpi::all_to_all(send.elem_nodes_id[f], recv.elem_nodes_id[f] );
      eckit::mpi::all_to_all(send.elem_part    [f], recv.elem_part    [f] );
    }
  }


public:
  Mesh& mesh;
  ArrayView<double,2> lonlat;
  ArrayView<gidx_t,1> glb_idx;
  ArrayView<int   ,1> part;
  IndexView<int   ,1> ridx;
  ArrayView<int   ,1> flags;
  std::vector< IndexView<int,   2> > elem_nodes;
  std::vector< ArrayView<int,   1> > elem_part;
  std::vector< ArrayView<gidx_t,1> > elem_glb_idx;
  std::vector<int> bdry_nodes;
  Node2Elem node_to_elem;
  Uid2Node uid2node;
  UniqueLonLat compute_uid;


public:
  BuildHaloHelper( Mesh& _mesh ): mesh(_mesh),
    compute_uid(mesh.nodes())
  {
    update();
  }

  void update()
  {
    compute_uid.update();
    Nodes& nodes         = mesh.nodes();
    lonlat   = ArrayView<double,2> ( nodes.lonlat() );
    glb_idx  = ArrayView<gidx_t,1> ( nodes.global_index() );
    part     = ArrayView<int   ,1> ( nodes.partition() );
    ridx     = IndexView<int   ,1> ( nodes.remote_index() );
    flags    = ArrayView<int   ,1> ( nodes.field("flags") );

    elem_nodes.  resize( mesh.nb_function_spaces() );
    elem_part.   resize( mesh.nb_function_spaces() );
    elem_glb_idx.resize( mesh.nb_function_spaces() );

    for(size_t f = 0; f < mesh.nb_function_spaces(); ++f)
    {
      FunctionSpace& elements = mesh.function_space(f);
      if( elements.metadata().get<long>("type") == Entity::ELEMS )
      {
        elem_nodes  [f] = IndexView<int,   2>( elements.field("nodes") );
        elem_part   [f] = ArrayView<int,   1>( elements.field("partition") );
        elem_glb_idx[f] = ArrayView<gidx_t,1>( elements.field("glb_idx") );
      }
    }
    compute_uid.update();
  }

  template< typename NodeContainer, typename ElementContainer >
  void fill_sendbuffer(Buffers& buf,const NodeContainer& nodes_uid, const ElementContainer& elems, const int p)
  {
    int nb_nodes = nodes_uid.size();
    buf.node_glb_idx[p].resize(nb_nodes);
    buf.node_part   [p].resize(nb_nodes);
    buf.node_ridx   [p].resize(nb_nodes);
    buf.node_flags  [p].resize(nb_nodes,Topology::NONE);
    buf.node_lonlat [p].resize(2*nb_nodes);

    int jnode=0;
    typename NodeContainer::iterator it;
    for( it=nodes_uid.begin(); it!=nodes_uid.end(); ++it, ++jnode )
    {
      uid_t uid = *it;

      Uid2Node::iterator found = uid2node.find( uid );
      if( found != uid2node.end() ) // Point exists inside domain
      {
        int node = found->second;
        buf.node_glb_idx[p][jnode]      = glb_idx(node);
        buf.node_part   [p][jnode]      = part   (node);
        buf.node_ridx   [p][jnode]      = ridx   (node);
        buf.node_lonlat [p][jnode*2+LON] = lonlat (node,LON);
        buf.node_lonlat [p][jnode*2+LAT] = lonlat (node,LAT);
        Topology::set(buf.node_flags[p][jnode],flags(node)|Topology::GHOST);
      }
      else
      {
        eckit::Log::warning() << "Node with uid " << uid << " needed by ["<<p<<"] was not found in ["<<eckit::mpi::rank()<<"]." << std::endl;
        ASSERT(false);
      }
    }

    for( int f=0; f<mesh.nb_function_spaces(); ++f )
    {
      FunctionSpace& elements = mesh.function_space(f);

      if( elements.metadata().get<long>("type") == Entity::ELEMS )
      {
        int nb_elem_nodes(elem_nodes[f].shape(1));
        buf.elem_glb_idx [f][p].resize( elems[f].size() );
        buf.elem_part    [f][p].resize( elems[f].size() );
        buf.elem_nodes_id[f][p].resize( elems[f].size()*nb_elem_nodes );

        ArrayView<uid_t,2> sfe_nodes_id_view( buf.elem_nodes_id[f][p].data(),
                                              make_shape(elems[f].size(),nb_elem_nodes).data() );

        for( int jelem=0; jelem<elems[f].size(); ++jelem )
        {
          int e = elems[f][jelem];
          buf.elem_part   [f][p][jelem] = elem_part[f][e];
          buf.elem_glb_idx[f][p][jelem] = compute_uid( elem_nodes[f][e] );

          for( int n=0; n<nb_elem_nodes; ++n)
            sfe_nodes_id_view(jelem,n) = compute_uid( elem_nodes[f](e,n) );
        }
      }
    }

  }

  template< typename NodeContainer, typename ElementContainer >
  void fill_sendbuffer(Buffers& buf,const NodeContainer& nodes_uid, const ElementContainer& elems, const PeriodicTransform& transform, int newflags, const int p)
  {
    int nb_nodes = nodes_uid.size();
    buf.node_glb_idx[p].resize(nb_nodes);
    buf.node_part   [p].resize(nb_nodes);
    buf.node_ridx   [p].resize(nb_nodes);
    buf.node_flags  [p].resize(nb_nodes,Topology::NONE);
    buf.node_lonlat [p].resize(2*nb_nodes);

    int jnode=0;
    typename NodeContainer::iterator it;
    for( it=nodes_uid.begin(); it!=nodes_uid.end(); ++it, ++jnode )
    {
      uid_t uid = *it;

      Uid2Node::iterator found = uid2node.find( uid );
      if( found != uid2node.end() ) // Point exists inside domain
      {
        int node = found->second;
        buf.node_part   [p][jnode]      = part   (node);
        buf.node_ridx   [p][jnode]      = ridx   (node);
        buf.node_lonlat [p][jnode*2+LON] = lonlat (node,LON);
        buf.node_lonlat [p][jnode*2+LAT] = lonlat (node,LAT);
        transform(&buf.node_lonlat[p][jnode*2],-1);
        // Global index of node is based on UID of destination
        buf.node_glb_idx[p][jnode]      = util::unique_lonlat(&buf.node_lonlat [p][jnode*2]);
        Topology::set(buf.node_flags[p][jnode],newflags);
      }
      else
      {
        eckit::Log::warning() << "Node with uid " << uid << " needed by ["<<p<<"] was not found in ["<<eckit::mpi::rank()<<"]." << std::endl;
        ASSERT(false);
      }
    }

    for( int f=0; f<mesh.nb_function_spaces(); ++f )
    {
      FunctionSpace& elements = mesh.function_space(f);

      if( elements.metadata().get<long>("type") == Entity::ELEMS )
      {
        int nb_elem_nodes(elem_nodes[f].shape(1));
        buf.elem_glb_idx [f][p].resize( elems[f].size() );
        buf.elem_part    [f][p].resize( elems[f].size() );
        buf.elem_nodes_id[f][p].resize( elems[f].size()*nb_elem_nodes );

        ArrayView<uid_t,2> sfe_nodes_id_view( buf.elem_nodes_id[f][p].data(),
                                              make_shape(elems[f].size(),nb_elem_nodes).data() );

        for( int jelem=0; jelem<elems[f].size(); ++jelem )
        {
          int e = elems[f][jelem];
          buf.elem_part   [f][p][jelem] = elem_part[f][e];

          std::vector<double> crds(nb_elem_nodes*2);
          for( int n=0; n<nb_elem_nodes; ++n)
          {
            double crd[] = { lonlat(elem_nodes[f](e,n),LON) , lonlat(elem_nodes[f](e,n),LAT) };
            transform(crd,-1);
            sfe_nodes_id_view(jelem,n) = util::unique_lonlat(crd);
            crds[n*2+LON] = crd[LON];
            crds[n*2+LAT] = crd[LAT];
          }
          // Global index of element is based on UID of destination
          buf.elem_glb_idx[f][p][jelem] = util::unique_lonlat( crds.data(), nb_elem_nodes );
        }
      }
    }
  }


  void add_nodes(Buffers& buf)
  {
    Nodes& nodes = mesh.nodes();
    int nb_nodes = nodes.size();

    // Nodes might be duplicated from different Tasks. We need to identify unique entries
    std::set<uid_t> node_uid;
    for( int jnode=0; jnode<nb_nodes; ++jnode )
    {
      node_uid.insert( compute_uid(jnode) );
    }
    std::vector< std::vector<int> > rfn_idx(eckit::mpi::size());
    for(size_t jpart = 0; jpart < eckit::mpi::size(); ++jpart)
    {
      rfn_idx[jpart].reserve(buf.node_glb_idx[jpart].size());
    }

    int nb_new_nodes=0;
    for(size_t jpart = 0; jpart < eckit::mpi::size(); ++jpart)
    {
      for(size_t n = 0; n < buf.node_glb_idx[jpart].size(); ++n)
      {
        double crd[] = { buf.node_lonlat[jpart][n*2+LON], buf.node_lonlat[jpart][n*2+LAT] };
        bool inserted = node_uid.insert( util::unique_lonlat(crd) ).second;
        if( inserted ) {
          rfn_idx[jpart].push_back(n);
        }
      }
      nb_new_nodes += rfn_idx[jpart].size();
    }

    // Resize nodes
    // ------------
    nodes.resize( nb_nodes+nb_new_nodes );
    flags   = ArrayView<int,   1>( nodes.field("flags") );
    glb_idx = ArrayView<gidx_t,1>( nodes.global_index() );
    part    = ArrayView<int,   1>( nodes.partition() );
    ridx    = IndexView<int,   1>( nodes.remote_index() );
    lonlat  = ArrayView<double,2>( nodes.lonlat() );

    compute_uid.update();

    // Add new nodes
    // -------------
    int new_node=0;
    for(size_t jpart = 0; jpart < eckit::mpi::size(); ++jpart)
    {
      for(size_t n = 0; n < rfn_idx[jpart].size(); ++n)
      {
        int loc_idx = nb_nodes+new_node;
        Topology::reset(flags(loc_idx),buf.node_flags[jpart][rfn_idx[jpart][n]]);
        glb_idx(loc_idx)    = buf.node_glb_idx [jpart][rfn_idx[jpart][n]];
        part   (loc_idx)    = buf.node_part    [jpart][rfn_idx[jpart][n]];
        ridx   (loc_idx)    = buf.node_ridx    [jpart][rfn_idx[jpart][n]];
        lonlat (loc_idx,LON) = buf.node_lonlat  [jpart][rfn_idx[jpart][n]*2+LON];
        lonlat (loc_idx,LAT) = buf.node_lonlat  [jpart][rfn_idx[jpart][n]*2+LAT];
        uid_t uid = compute_uid(loc_idx);

        // make sure new node was not already there
        Uid2Node::iterator found = uid2node.find(uid);
        if( found != uid2node.end() )
        {
          int other = found->second;
          std::stringstream msg;
          msg << "New node with uid " << uid << ":\n"  << glb_idx(loc_idx)
              << "("<<lonlat(loc_idx,LON)<<","<<lonlat(loc_idx,LAT)<<")\n";
          msg << "Existing already loc "<< other << "  :  " << glb_idx(other)
              << "("<<lonlat(other,LON)<<","<<lonlat(other,LAT)<<")\n";
          throw eckit::SeriousBug(msg.str(),Here());
        }
        uid2node[ uid ] = nb_nodes+new_node;
        ++new_node;
      }
    }
  }


  void add_elements(Buffers& buf)
  {
    for(size_t f = 0; f < mesh.nb_function_spaces(); ++f)
    {
      FunctionSpace& elements = mesh.function_space(f);
      if( elements.metadata().get<long>("type") == Entity::ELEMS )
      {
        // Elements might be duplicated from different Tasks. We need to identify unique entries
        std::set<uid_t> elem_uid;
        int nb_elems = elements.shape(0);
        for( int jelem=0; jelem<nb_elems; ++jelem )
        {
          elem_uid.insert( compute_uid(elem_nodes[f][jelem]) );
        }

        std::vector< std::vector<int> > received_new_elems(eckit::mpi::size());
        for(size_t jpart = 0; jpart < eckit::mpi::size(); ++jpart)
        {
          received_new_elems[jpart].reserve(buf.elem_glb_idx[f][jpart].size());
        }

        int nb_new_elems=0;
        for(size_t jpart = 0; jpart < eckit::mpi::size(); ++jpart)
        {
          for(size_t e = 0; e < buf.elem_glb_idx[f][jpart].size(); ++e)
          {
            bool inserted = elem_uid.insert( buf.elem_glb_idx[f][jpart][e] ).second;
            if( inserted )
            {
              received_new_elems[jpart].push_back(e);
            }
          }
          nb_new_elems += received_new_elems[jpart].size();
        }

        // Resize elements
        // ---------------
        elements.resize( make_shape( nb_elems+nb_new_elems, FunctionSpace::UNDEF_VARS ) );
        elem_glb_idx[f] = ArrayView<gidx_t,1>( elements.field("glb_idx") );
        elem_nodes[f]   = IndexView<int,2>( elements.field("nodes")   );
        elem_part[f]    = ArrayView<int,1>( elements.field("partition")   );

        // Add new elems
        // -------------
        int nb_nodes_per_elem = elem_nodes[f].shape(1);
        int new_elem=0;
        for(size_t jpart = 0; jpart < eckit::mpi::size(); ++jpart)
        {
          for(size_t e = 0; e < received_new_elems[jpart].size(); ++e)
          {
            int jelem = received_new_elems[jpart][e];
            elem_glb_idx[f](nb_elems+new_elem)   = buf.elem_glb_idx[f][jpart][jelem];
            elem_part   [f](nb_elems+new_elem)   = buf.elem_part[f][jpart][jelem];
            for( int n=0; n<nb_nodes_per_elem; ++n )
              elem_nodes[f](nb_elems+new_elem,n) = uid2node[ buf.elem_nodes_id[f][jpart][jelem*nb_nodes_per_elem+n] ];
            ++new_elem;
          }
        }
      }
    }
  }

  void add_buffers(Buffers& buf)
  {
    add_nodes(buf);
    add_elements(buf);
    update();
  }

};


void increase_halo_interior( BuildHaloHelper& helper )
{
  if (helper.node_to_elem.size() == 0 )
    build_lookup_node2elem(helper.mesh,helper.node_to_elem);

  if( helper.uid2node.size() == 0 )
    build_lookup_uid2node(helper.mesh,helper.uid2node);


  // All buffers needed to move elements and nodes
  BuildHaloHelper::Buffers sendmesh(helper.mesh);
  BuildHaloHelper::Buffers recvmesh(helper.mesh);

  // 1) Find boundary nodes of this partition:

  accumulate_partition_bdry_nodes(helper.mesh,helper.bdry_nodes);
  const std::vector<int>& bdry_nodes = helper.bdry_nodes;

  // 2) Communicate uid of these boundary nodes to other partitions

  std::vector<uid_t> send_bdry_nodes_uid(bdry_nodes.size());
  for(size_t jnode = 0; jnode < bdry_nodes.size(); ++jnode)
    send_bdry_nodes_uid[jnode] = helper.compute_uid(bdry_nodes[jnode]);

  atlas::mpi::Buffer<uid_t,1> recv_bdry_nodes_uid_from_parts;
  eckit::mpi::all_gather(send_bdry_nodes_uid,recv_bdry_nodes_uid_from_parts);

  for (size_t jpart = 0; jpart < eckit::mpi::size(); ++jpart)
  {
    // 3) Find elements and nodes completing these elements in
    //    other tasks that have my nodes through its UID

    ArrayView<uid_t,1> recv_bdry_nodes_uid = recv_bdry_nodes_uid_from_parts[jpart];

    std::vector< std::vector<int> > found_bdry_elems;
    std::set< uid_t > found_bdry_nodes_uid;

    accumulate_elements(helper.mesh,recv_bdry_nodes_uid,
                        helper.uid2node,
                        helper.node_to_elem,
                        found_bdry_elems,
                        found_bdry_nodes_uid);

    // 4) Fill node and element buffers to send back
    helper.fill_sendbuffer(sendmesh,found_bdry_nodes_uid,found_bdry_elems,jpart);
  }

  // 5) Now communicate all buffers
  helper.all_to_all(sendmesh,recvmesh);

  // 6) Adapt mesh
  helper.add_buffers(recvmesh);
}

class PeriodicPoints {
public:
  PeriodicPoints(Mesh& mesh, int flag, size_t N)
  {
    flag_ = flag;
    N_ = N;
    flags_ = ArrayView<int,1> ( mesh.nodes().field("flags") );
  }

  bool operator()(int j) const
  {
    if( j>=N_ ) return false;
    if( Topology::check(flags_(j),flag_) ) return true;
    return false;
  }
private:
  int N_;
  int flag_;
  ArrayView<int,1> flags_;
};

void increase_halo_periodic( BuildHaloHelper& helper, const PeriodicPoints& periodic_points, const PeriodicTransform& transform, int newflags)
{
  if (helper.node_to_elem.size() == 0 )
    build_lookup_node2elem(helper.mesh,helper.node_to_elem);

  if( helper.uid2node.size() == 0 )
    build_lookup_uid2node(helper.mesh,helper.uid2node);

  // All buffers needed to move elements and nodes
  BuildHaloHelper::Buffers sendmesh(helper.mesh);
  BuildHaloHelper::Buffers recvmesh(helper.mesh);

  // 1) Find boundary nodes of this partition:

  if( ! helper.bdry_nodes.size() )
    accumulate_partition_bdry_nodes(helper.mesh,helper.bdry_nodes);

  std::vector<int> bdry_nodes = filter_nodes(helper.bdry_nodes,periodic_points);

  //DEBUG_VAR(bdry_nodes.size());

  // 2) Compute transformed uid of these boundary nodes and send to other partitions

  std::vector<uid_t> send_bdry_nodes_uid(bdry_nodes.size());
  for(size_t jnode = 0; jnode < bdry_nodes.size(); ++jnode)
  {
    double crd[] = { helper.lonlat(bdry_nodes[jnode],LON), helper.lonlat(bdry_nodes[jnode],LAT) };
    transform(crd,+1);
    send_bdry_nodes_uid[jnode] = util::unique_lonlat(crd);
  }
  atlas::mpi::Buffer<uid_t,1> recv_bdry_nodes_uid_from_parts;
  eckit::mpi::all_gather(send_bdry_nodes_uid,recv_bdry_nodes_uid_from_parts);

  for (size_t jpart = 0; jpart < eckit::mpi::size(); ++jpart)
  {
    // 3) Find elements and nodes completing these elements in
    //    other tasks that have my nodes through its UID

    ArrayView<uid_t,1> recv_bdry_nodes_uid = recv_bdry_nodes_uid_from_parts[jpart];

    std::vector< std::vector<int> > found_bdry_elems;
    std::set< uid_t > found_bdry_nodes_uid;

    accumulate_elements(helper.mesh,recv_bdry_nodes_uid,
                        helper.uid2node,
                        helper.node_to_elem,
                        found_bdry_elems,
                        found_bdry_nodes_uid);

    //DEBUG("elements from part["<<jpart<<"] = " << found_bdry_elems.size());

    // 4) Fill node and element buffers to send back
    helper.fill_sendbuffer(sendmesh,found_bdry_nodes_uid,found_bdry_elems,transform,newflags,jpart);
  }

  // 5) Now communicate all buffers
  helper.all_to_all(sendmesh,recvmesh);

  // 6) Adapt mesh

  helper.add_buffers(recvmesh);

}




void build_halo(Mesh& mesh, int nb_elems )
{
  int jhalo = 0;
  mesh.metadata().get("halo",jhalo);

  for( ; jhalo<nb_elems; ++jhalo )
  {
    size_t nb_nodes_before_halo_increase = mesh.nodes().size();

    BuildHaloHelper helper(mesh);

    increase_halo_interior( helper );

    PeriodicPoints westpts(mesh,Topology::PERIODIC|Topology::WEST,nb_nodes_before_halo_increase);

    increase_halo_periodic( helper, westpts, WestEast(), Topology::PERIODIC|Topology::WEST|Topology::GHOST );

    PeriodicPoints eastpts(mesh,Topology::PERIODIC|Topology::EAST,nb_nodes_before_halo_increase);

    increase_halo_periodic( helper, eastpts, EastWest(), Topology::PERIODIC|Topology::EAST|Topology::GHOST );

    std::stringstream ss;
    ss << "nb_nodes_including_halo["<<jhalo+1<<"]";
    mesh.metadata().set(ss.str(),mesh.nodes().size());
  }

  mesh.metadata().set("halo",nb_elems);
}








// ------------------------------------------------------------------
// C wrapper interfaces to C++ routines

void atlas__build_halo ( Mesh* mesh, int nb_elems ) {
#undef ATLAS_ERROR_HANDLING
#define ATLAS_ERROR_HANDLING(x) x
  ATLAS_ERROR_HANDLING( build_halo(*mesh, nb_elems) );
}

// ------------------------------------------------------------------

} // namespace actions
} // namespace atlas

