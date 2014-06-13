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
#include "atlas/mesh/Mesh.hpp"
#include "atlas/mesh/FunctionSpace.hpp"
#include "atlas/mesh/Field.hpp"
#include "atlas/actions/BuildHalo.hpp"
#include "atlas/mesh/Parameters.hpp"
#include "atlas/mesh/ArrayView.hpp"
#include "atlas/mesh/IndexView.hpp"
#include "atlas/mesh/Array.hpp"
#include "atlas/mesh/Util.hpp"

namespace atlas {
namespace actions {

namespace buildhalo_detail {


int EAST = BC::EAST;
int WEST = BC::WEST;

}

using namespace buildhalo_detail;

// ------------------------------------------------------------------

void increase_halo( Mesh& mesh );

void build_halo(Mesh& mesh, int nb_elems )
{
  for( int jhalo=0; jhalo<nb_elems; ++jhalo )
  {
    increase_halo( mesh );
  }
}

void increase_halo( Mesh& mesh )
{
  FunctionSpace& nodes         = mesh.function_space( "nodes" );
  ArrayView<double,2> latlon   ( nodes.field( "coordinates"    ) );
  ArrayView<int   ,1> glb_idx  ( nodes.field( "glb_idx"        ) );
  ArrayView<int   ,1> part     ( nodes.field( "partition"      ) );
  IndexView<int   ,1> ridx     ( nodes.field( "remote_idx"     ) );

  int nb_nodes = nodes.extents()[0];

  std::vector< std::vector< ElementRef > > node_to_elem(nb_nodes);
  
  std::vector< IndexView<int,2> > elem_nodes( mesh.nb_function_spaces() );
  std::vector< ArrayView<int,1> > elem_part ( mesh.nb_function_spaces() );
  std::vector< ArrayView<int,1> > elem_glb_idx ( mesh.nb_function_spaces() );

  for( int func_space_idx=0; func_space_idx<mesh.nb_function_spaces(); ++func_space_idx)
  {
    FunctionSpace& elements = mesh.function_space(func_space_idx);
    if( elements.metadata<int>("type") == Entity::ELEMS )
    {
      elem_nodes  [func_space_idx] = IndexView<int,2>( elements.field("nodes") );
      elem_part   [func_space_idx] = ArrayView<int,1>( elements.field("partition") );
      elem_glb_idx[func_space_idx] = ArrayView<int,1>( elements.field("glb_idx") );
      int nb_elems = elem_nodes[func_space_idx].extents()[0];
      int nb_nodes_per_elem = elem_nodes[func_space_idx].extents()[1];
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


  int extents[] = {nb_faces,2};
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

  std::map<int,int> node_uid_to_loc;
  for( int jnode=0; jnode<nb_nodes; ++jnode )
  {
    LatLonPoint ll(latlon[jnode]);
    if( node_uid_to_loc.count(ll.uid()) > 0 )
    {
      int other = node_uid_to_loc[ll.uid()];
      std::stringstream msg;
      msg << "Node uid: " << ll.uid() << "   " << glb_idx(jnode)
          << " (" << latlon(jnode,XX) <<","<< latlon(jnode,YY)<<")  has already been added as node "
          << glb_idx(other) << " (" << latlon(other,XX) <<","<< latlon(other,YY)<<")";
      throw eckit::SeriousBug(msg.str(),Here());
    }
    node_uid_to_loc[ll.uid()] = jnode;
  }

  int nb_bdry_nodes = bdry_nodes.size();
  Array<int> arr_bdry_nodes_id(nb_bdry_nodes,3);
  ArrayView<int,2> bdry_nodes_id(arr_bdry_nodes_id);
  ASSERT( bdry_nodes_id.extents()[0] == nb_bdry_nodes );
  ASSERT( bdry_nodes_id.extents()[1] == 3);

  for( int jnode=0; jnode<nb_bdry_nodes; ++jnode )
  {
    LatLonPoint ll( latlon[bdry_nodes[jnode]] );
    bdry_nodes_id(jnode,0) = ll.x;
    bdry_nodes_id(jnode,1) = ll.y;
    bdry_nodes_id(jnode,2) = glb_idx( bdry_nodes[jnode] );
  }



  std::vector<int> recvcounts( MPL::size() );
  std::vector<int> recvdispls( MPL::size() );
  int sendcnt = bdry_nodes_id.total_size();
  ASSERT( sendcnt == nb_bdry_nodes*3 );
  MPL_CHECK_RESULT( MPI_Allgather( &sendcnt,          1, MPI_INT,
                                   recvcounts.data(), 1, MPI_INT, MPI_COMM_WORLD ) );

  recvdispls[0] = 0;
  int recvcnt = recvcounts[0];
  for( int jpart=1; jpart<MPL::size(); ++jpart )
  {
    recvdispls[jpart] = recvdispls[jpart-1] + recvcounts[jpart-1];
    recvcnt += recvcounts[jpart];
  }
  std::vector<int> recvbuf(recvcnt);

  MPL_CHECK_RESULT( MPI_Allgatherv( bdry_nodes_id.data(), sendcnt, MPI_INT,
                    recvbuf.data(), recvcounts.data(), recvdispls.data(),
                    MPI_INT, MPI_COMM_WORLD) );


  // sfn stands for "send_found_nodes"
  std::vector< std::vector<int>    > sfn_part( MPL::size() );
  std::vector< std::vector<int>    > sfn_ridx( MPL::size() );
  std::vector< std::vector<int>    > sfn_glb_idx( MPL::size() );
  std::vector< std::vector<double> > sfn_latlon ( MPL::size() );
  // sfn stands for "send_found_elems"
  std::vector< std::vector< std::vector<int> > >
      sfe_glb_idx ( mesh.nb_function_spaces(), std::vector< std::vector<int> >( MPL::size() ) );
  std::vector< std::vector< std::vector<int> > >
      sfe_nodes_id( mesh.nb_function_spaces(), std::vector< std::vector<int> >( MPL::size() ) );
  std::vector< std::vector< std::vector<int> > >
      sfe_part    ( mesh.nb_function_spaces(), std::vector< std::vector<int> >( MPL::size() ) );

  for (int jpart=0; jpart<MPL::size(); ++jpart)
  {
    ArrayView<int,2> recv_bdry_nodes_id( recvbuf.data()+recvdispls[jpart],
                                         Extents( recvcounts[jpart]/3, 3 ).data() );
    int recv_nb_bdry_nodes = recv_bdry_nodes_id.extents()[0];

    // Find elements that have these nodes
    // In order to do this, check the node_to_elem list
    // Warning: only add elements that are owned!


    std::vector< std::set< std::pair<int,int> > > found_bdry_elements_set( mesh.nb_function_spaces() );
    for( int jrecv=0; jrecv<recv_nb_bdry_nodes; ++jrecv )
    {
      ASSERT( recv_bdry_nodes_id.extents()[1] == 3 );
      int recv_x       = recv_bdry_nodes_id(jrecv,0);
      int recv_y       = recv_bdry_nodes_id(jrecv,1);
      int recv_glb_idx = recv_bdry_nodes_id(jrecv,2);

      // Only search for nodes in the interior latlon domain
      LatLonPoint ll(recv_x,recv_y);
      if     ( ll.x <= WEST) { while(ll.x <= WEST) { ll.x += EAST; } }
      else if( ll.x >= EAST) { while(ll.x >= EAST) { ll.x -= EAST; } }

      int recv_uid = ll.uid();
      int loc = -1;
      std::map<int,int>::iterator found = node_uid_to_loc.find(recv_uid);
      if( found != node_uid_to_loc.end() )
      {
        loc = found->second;
        if( MPL::rank() == jpart && glb_idx(loc) == recv_glb_idx )
        {
          loc = -1;
        }
      }

      if( loc != -1 )
      {
        for( int jelem=0; jelem<node_to_elem[loc].size(); ++jelem )
        {
          int f = node_to_elem[loc][jelem].f;
          int e = node_to_elem[loc][jelem].e;

          if( elem_part[f](e) == MPL::rank() )
          {
            found_bdry_elements_set[f].insert( std::make_pair( e, recv_x - ll.x ) );
          }
        }
      }
    }
    std::vector< std::vector<int> > found_bdry_elements( mesh.nb_function_spaces() );
    std::vector< std::vector<int> > found_bdry_elements_coord_transform( mesh.nb_function_spaces() );

    std::vector<int> nb_found_bdry_elems( mesh.nb_function_spaces(), 0 );

    for( int f=0; f<mesh.nb_function_spaces(); ++f )
    {
      nb_found_bdry_elems[f] = found_bdry_elements_set[f].size();

      found_bdry_elements[f].resize(nb_found_bdry_elems[f]);
      found_bdry_elements_coord_transform[f].resize(nb_found_bdry_elems[f]);
      int jelem=0;
      std::set< std::pair<int,int> >::iterator it=found_bdry_elements_set[f].begin();
      for( ; it!=found_bdry_elements_set[f].end(); ++it, ++jelem )
      {
        found_bdry_elements[f][jelem] = it->first;
        found_bdry_elements_coord_transform[f][jelem] = it->second;
      }
      nb_found_bdry_elems[f] = found_bdry_elements[f].size();
    }

    // Collect all nodes needed to complete the element
    std::set<LatLonPoint> found_bdry_nodes_id_set;
    {
      for( int f=0; f<mesh.nb_function_spaces(); ++f )
      {
        for( int jelem=0; jelem<nb_found_bdry_elems[f]; ++jelem )
        {
          int e = found_bdry_elements[f][jelem];
          int nb_elem_nodes = elem_nodes[f].extents()[1];
          for( int n=0; n<nb_elem_nodes; ++n )
          {
            int x = microdeg( latlon( elem_nodes[f](e,n), XX) ) + found_bdry_elements_coord_transform[f][jelem];
            int y = microdeg( latlon( elem_nodes[f](e,n), YY) );
            found_bdry_nodes_id_set.insert( LatLonPoint(x,y) );
          }
        }
      }

      // Remove nodes we already have received, as we won't need to send them
      for( int jrecv=0; jrecv<recv_nb_bdry_nodes; ++jrecv)
      {
        int x = recv_bdry_nodes_id(jrecv,0);
        int y = recv_bdry_nodes_id(jrecv,1);
        found_bdry_nodes_id_set.erase( LatLonPoint(x,y) ) ;
      }
    }
    int nb_found_bdry_nodes = found_bdry_nodes_id_set.size();

    sfn_glb_idx[jpart].resize(nb_found_bdry_nodes);
    sfn_part[jpart].resize(nb_found_bdry_nodes);
    sfn_ridx[jpart].resize(nb_found_bdry_nodes);
    sfn_latlon[jpart].resize(2*nb_found_bdry_nodes);

    // Fill buffers to send
    {
      int jnode=0;
      std::set<LatLonPoint>::iterator it=found_bdry_nodes_id_set.begin();
      for( ; it!=found_bdry_nodes_id_set.end(); ++it, ++jnode )
      {
        int uid = it->uid();
        std::map<int,int>::iterator found = node_uid_to_loc.find( uid );
        if( found != node_uid_to_loc.end() ) // Point exists inside domain
        {
          int loc = found->second;
          sfn_glb_idx[jpart][jnode]      = glb_idx(loc);
          sfn_part   [jpart][jnode]      = part   (loc);
          sfn_ridx   [jpart][jnode]      = ridx   (loc);
          sfn_latlon [jpart][jnode*2+XX] = latlon (loc,XX);
          sfn_latlon [jpart][jnode*2+YY] = latlon (loc,YY);
        }
        else // periodic point
        {
          LatLonPoint ll(it->x,it->y);
          // Here ll has the new coords
          if     ( ll.x <= WEST) { while(ll.x <= WEST) { ll.x += EAST; } }
          else if( ll.x >= EAST) { while(ll.x >= EAST) { ll.x -= EAST; } }
          // Now ll has the coords of periodic point
          int pid = ll.uid();
          found = node_uid_to_loc.find( pid );
          ASSERT( found != node_uid_to_loc.end() );
          int loc = found->second;
          sfn_glb_idx[jpart][jnode]      = uid;
          sfn_part   [jpart][jnode]      = part(loc);
          sfn_ridx   [jpart][jnode]      = ridx(loc);
          sfn_latlon [jpart][jnode*2+XX] = latlon(loc,XX) + (it->x - ll.x)/EAST * 2.*M_PI;
          sfn_latlon [jpart][jnode*2+YY] = latlon(loc,YY);
        }
      }
    }

    for( int f=0; f<mesh.nb_function_spaces(); ++f )
    {
      int nb_elem_nodes = elem_nodes[f].extents()[1];
      sfe_glb_idx [f][jpart].resize( nb_found_bdry_elems[f] );
      sfe_part    [f][jpart].resize( nb_found_bdry_elems[f] );
      sfe_nodes_id[f][jpart].resize( nb_found_bdry_elems[f]*nb_elem_nodes );
      ArrayView<int,2> sfe_nodes_id_view( sfe_nodes_id[f][jpart].data(),
                                          Extents(nb_found_bdry_elems[f],nb_elem_nodes).data() );
      for( int jelem=0; jelem<nb_found_bdry_elems[f]; ++jelem )
      {
        int e = found_bdry_elements[f][jelem];
        int xper = found_bdry_elements_coord_transform[f][jelem];
        sfe_part[f][jpart][jelem]    = elem_part[f][e];
        double centroid[2];
        centroid[XX] = 0.;
        centroid[YY] = 0.;
        for( int n=0; n<nb_elem_nodes; ++n)
        {
          double x, y;
          x = latlon(elem_nodes[f](e,n),XX) + xper/EAST*2.*M_PI;
          y = latlon(elem_nodes[f](e,n),YY);
          centroid[XX] += x;
          centroid[YY] += y;
          sfe_nodes_id_view(jelem,n) = LatLonPoint( x, y ).uid();
        }
        centroid[XX] /= static_cast<double>(nb_elem_nodes);
        centroid[YY] /= static_cast<double>(nb_elem_nodes);
        sfe_glb_idx[f][jpart][jelem] = LatLonPoint( centroid[XX], centroid[YY ]).uid() ;
      }
    }
  }

  // Now communicate all found fields back

  //    rfn stands for "recv_found_nodes"
  std::vector< std::vector<int> >    rfn_glb_idx(MPL::size());
  std::vector< std::vector<int> >    rfn_part(MPL::size());
  std::vector< std::vector<int>    > rfn_ridx( MPL::size() );
  std::vector< std::vector<double> > rfn_latlon(MPL::size());
  //    rfe stands for "recv_found_elems"
  std::vector< std::vector< std::vector<int> > >
      rfe_glb_idx ( mesh.nb_function_spaces(), std::vector< std::vector<int> >( MPL::size() ) );
  std::vector< std::vector< std::vector<int> > >
      rfe_part    ( mesh.nb_function_spaces(), std::vector< std::vector<int> >( MPL::size() ) );
  std::vector< std::vector< std::vector<int> > >
      rfe_nodes_id( mesh.nb_function_spaces(), std::vector< std::vector<int> >( MPL::size() ) );

  MPL::Alltoall(sfn_glb_idx,  rfn_glb_idx);
  MPL::Alltoall(sfn_part,     rfn_part);
  MPL::Alltoall(sfn_ridx,     rfn_ridx);
  MPL::Alltoall(sfn_latlon,   rfn_latlon);
  for( int f=0; f<mesh.nb_function_spaces(); ++f )
  {
    MPL::Alltoall(sfe_glb_idx [f], rfe_glb_idx [f] );
    MPL::Alltoall(sfe_nodes_id[f], rfe_nodes_id[f] );
    MPL::Alltoall(sfe_part    [f], rfe_part    [f] );
  }

  // We now have everything we need in rfe_ and rfn_ vectors
  // Now adapt the mesh

  // Nodes might be duplicated from different Tasks. We need to identify unique entries
  std::set<int> node_uid;
  for( int jnode=0; jnode<nb_nodes; ++jnode )
  {
    node_uid.insert( LatLonPoint( latlon[jnode] ).uid() );
  }
  std::vector< std::vector<int> > rfn_idx(MPL::size());
  for( int jpart=0; jpart<MPL::size(); ++jpart )
  {
    rfn_idx[jpart].reserve(rfn_glb_idx[jpart].size());
  }

  int nb_new_nodes=0;
  for( int jpart=0; jpart<MPL::size(); ++jpart )
  {
    for( int n=0; n<rfn_glb_idx[jpart].size(); ++n )
    {
      double x = rfn_latlon[jpart][n*2+XX];
      double y = rfn_latlon[jpart][n*2+YY];
      bool inserted = node_uid.insert( LatLonPoint( x, y ).uid() ).second;
      if( inserted )
        rfn_idx[jpart].push_back(n);
    }
    nb_new_nodes += rfn_idx[jpart].size();
  }


  nodes.resize( Extents( nb_nodes+nb_new_nodes, Field::UNDEF_VARS ) );
  glb_idx = ArrayView<int,   1>( nodes.field("glb_idx") );
  part    = ArrayView<int,   1>( nodes.field("partition") );
  ridx    = IndexView<int,   1>( nodes.field("remote_idx") );
  latlon  = ArrayView<double,2>( nodes.field("coordinates") );


  int new_node=0;
  for( int jpart=0; jpart<MPL::size(); ++jpart )
  {
    for( int n=0; n<rfn_idx[jpart].size(); ++n )
    {
      int loc_idx = nb_nodes+new_node;

      glb_idx(loc_idx)    = rfn_glb_idx [jpart][rfn_idx[jpart][n]];
      part   (loc_idx)    = rfn_part    [jpart][rfn_idx[jpart][n]];
      ridx   (loc_idx)    = rfn_ridx    [jpart][rfn_idx[jpart][n]];
      latlon (loc_idx,XX) = rfn_latlon  [jpart][rfn_idx[jpart][n]*2+XX];
      latlon (loc_idx,YY) = rfn_latlon  [jpart][rfn_idx[jpart][n]*2+YY];
      int uid = LatLonPoint( latlon[loc_idx] ).uid();

      // make sure new node was not already there
      std::map<int,int>::iterator found = node_uid_to_loc.find(uid);
      if( found != node_uid_to_loc.end() )
      {
        int other = found->second;
        std::stringstream msg;
        msg << "New node with uid " << uid << ":\n"  << glb_idx(loc_idx)
            << "("<<latlon(loc_idx,XX)<<","<<latlon(loc_idx,YY)<<")\n";
        msg << "Existing already loc "<< other << "  :  " << glb_idx(other)
            << "("<<latlon(other,XX)<<","<<latlon(other,YY)<<")\n";
        throw eckit::SeriousBug(msg.str(),Here());
      }
      node_uid_to_loc[ uid ] = nb_nodes+new_node;
      ++new_node;
    }
  }

  for( int f=0; f<mesh.nb_function_spaces(); ++f )
  {
    FunctionSpace& elements = mesh.function_space(f);
    if( elements.metadata<int>("type") == Entity::ELEMS )
    {
      ComputeUniqueElementIndex compute_uid( nodes );

      std::set<int> elem_uid;
      int nb_elems = elements.extents()[0];
      for( int jelem=0; jelem<nb_elems; ++jelem )
      {
        elem_uid.insert( compute_uid(elem_nodes[f][jelem]) );
      }

      std::vector< std::vector<int> > rfe_unique_idx(MPL::size());
      for( int jpart=0; jpart<MPL::size(); ++jpart )
      {
        rfe_unique_idx[jpart].reserve(rfe_glb_idx[f][jpart].size());
      }

      int nb_new_elems=0;
      for( int jpart=0; jpart<MPL::size(); ++jpart )
      {
        for( int e=0; e<rfe_glb_idx[f][jpart].size(); ++e )
        {
          bool inserted = elem_uid.insert( rfe_glb_idx[f][jpart][e] ).second;
          if( inserted )
          {
            rfe_unique_idx[jpart].push_back(e);
          }
        }
        nb_new_elems += rfe_unique_idx[jpart].size();
      }

      int nb_nodes_per_elem = elem_nodes[f].extents()[1];
      elements.resize( Extents( nb_elems+nb_new_elems, Field::UNDEF_VARS ) );
      elem_glb_idx[f] = ArrayView<int,1>( elements.field("glb_idx") );
      elem_nodes[f]   = IndexView<int,2>( elements.field("nodes")   );
      elem_part[f]    = ArrayView<int,1>( elements.field("partition")   );
      int new_elem=0;
      for( int jpart=0; jpart<MPL::size(); ++jpart )
      {
        for( int e=0; e<rfe_unique_idx[jpart].size(); ++e )
        {
          int jelem = rfe_unique_idx[jpart][e];
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

// ------------------------------------------------------------------
// C wrapper interfaces to C++ routines

void atlas__build_halo ( Mesh* mesh, int nb_elems ) {
  build_halo(*mesh, nb_elems);
}

// ------------------------------------------------------------------

} // namespace actions
} // namespace atlas

