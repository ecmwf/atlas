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
#include "atlas/mesh/ArrayView.hpp"
#include "atlas/mesh/IndexView.hpp"
#include "atlas/mesh/Array.hpp"

namespace atlas {
namespace actions {

ArrayView<int,1> build_nodes_glb_idx_serial( FunctionSpace& nodes );
ArrayView<int,1> build_nodes_partition( FunctionSpace& nodes );
IndexView<int,1> build_nodes_remote_loc_idx( FunctionSpace& nodes );

int microdeg( const double& deg )
{
  return static_cast<int>(deg*1.e6);
}

// Node struct that holds the longitude and latitude in millidegrees (integers)
// This structure is used in sorting algorithms, and uses less memory than
// if x and y were in double precision.
struct LatLon
{
  LatLon() {}
  LatLon( int x_, int y_ )
  {
    x = x_;
    y = y_;
  }
  LatLon( const ArrayView<int,1>& coord )
  {
    x = coord[XX];
    y = coord[YY];
  }
  LatLon( const ArrayView<double,1>& coord )
  {
    x = microdeg(coord[XX]);
    y = microdeg(coord[YY]);
  }
  int x, y;
  bool operator < (const LatLon& other) const
  {
    if( y > other.y  ) return true;
    if( y == other.y ) return (x < other.x);
    return false;
  }
};

// ------------------------------------------------------------------

void build_parallel_fields( Mesh& mesh )
{
  ASSERT( mesh.has_function_space("nodes") );

  FunctionSpace& nodes = mesh.function_space("nodes");

  ArrayView<int,1> nodes_glb_idx;
  if( nodes.has_field("glb_idx") )
    nodes_glb_idx = ArrayView<int,1> ( nodes.field("glb_idx") );
  else if( MPL::size() > 1 )
    throw eckit::SeriousBug( "parallel mesh needs nodes glb_idx", Here() );
  else
    nodes_glb_idx = build_nodes_glb_idx_serial( nodes );

  ArrayView<int,1> nodes_part;
  if( nodes.has_field("partition") )
    nodes_part = ArrayView<int,1> ( nodes.field("partition") );
  else
    nodes_part = build_nodes_partition( nodes );

  IndexView<int,1> nodes_loc_idx;
  if( nodes.has_field("remote_loc_idx") )
    nodes_loc_idx = IndexView<int,1> ( nodes.field("remote_loc_idx") );
  else
    nodes_loc_idx = build_nodes_remote_loc_idx( nodes );
}

// ------------------------------------------------------------------

ArrayView<int,1> build_nodes_glb_idx_serial( FunctionSpace& nodes )
{
  ArrayView<int,1> glb_idx ( nodes.create_field<int>("glb_idx",1) );
  for( int jnode=0; jnode<glb_idx.extents()[0]; ++jnode )
    glb_idx[jnode] = jnode+1;
  return glb_idx;
}

// ------------------------------------------------------------------

IndexView<int,1> build_nodes_remote_loc_idx( FunctionSpace& nodes )
{
  int mypart = MPL::rank();
  IndexView<int,   1> loc_idx ( nodes.create_field<int>("remote_loc_idx",1) );
  ArrayView<int,   1> part    ( nodes.field("partition")   );
  ArrayView<double,2> latlon  ( nodes.field("coordinates") );

  std::vector< std::vector<int> > send_needed( MPL::size() );
  std::vector< std::vector<int> > recv_needed( MPL::size() );

  int nb_nodes = nodes.extents()[0];
  int sendcnt=0;
  std::map<LatLon,int> lookup;
  for( int jnode=0; jnode<nb_nodes; ++jnode )
  {
    LatLon ll(latlon[jnode]);
    if( part(jnode)==mypart )
    {
      lookup[ ll ] = jnode;
      loc_idx(jnode) = jnode;
    }
    else
    {
      send_needed[part(jnode)].push_back( ll.x  );
      send_needed[part(jnode)].push_back( ll.y  );
      send_needed[part(jnode)].push_back( jnode );
      sendcnt++;
    }
  }

  MPL::Alltoall( send_needed, recv_needed );

  std::vector< std::vector<int> > send_found( MPL::size() );
  std::vector< std::vector<int> > recv_found( MPL::size() );

  for( int jpart=0; jpart<MPL::size(); ++jpart )
  {
    ArrayView<int,2> recv_node( recv_needed[jpart].data(), Extents(recv_needed[jpart].size()/3,3).data() );
    for( int jnode=0; jnode<recv_node.extents()[0]; ++jnode )
    {
      LatLon ll( recv_node[jnode] );
      if( lookup.count(ll) )
      {
        send_found[jpart].push_back( recv_node(jnode,2) );
        send_found[jpart].push_back( lookup[ll] );
      }
    }
  }

  MPL::Alltoall( send_found, recv_found );

  for( int jpart=0; jpart<MPL::size(); ++jpart )
  {
    ArrayView<int,2> recv_node( recv_found[jpart].data(), Extents(recv_found[jpart].size()/2,2).data() );
    for( int jnode=0; jnode<recv_node.extents()[0]; ++jnode )
    {
      loc_idx( recv_node(jnode,0) ) = recv_node(jnode,1);
    }
  }

  return loc_idx;
}

// ------------------------------------------------------------------

ArrayView<int,1> build_nodes_partition( FunctionSpace& nodes )
{
  ArrayView<int,1> part ( nodes.create_field<int>("partition",1) );
  if( MPL::size() == 1 )
    part = MPL::rank();
  else
  {
    if( nodes.has_field("proc") )
    {
      ArrayView<int,1> proc( nodes.field("proc") );
      for( int jnode=0; jnode<nodes.extents()[0]; ++jnode )
        part(jnode) = proc(jnode);
    }
    else
      NOTIMP;
  }
  return part;
}

// ------------------------------------------------------------------

void make_periodic( Mesh& mesh )
{
  int mypart = MPL::rank();

  FunctionSpace& nodes = mesh.function_space("nodes");

  IndexView<int,1> loc_idx ( nodes.field("remote_loc_idx") );
  ArrayView<int,1> part    ( nodes.field("partition")      );
  //ArrayView<int,1> glb_idx ( nodes.field("glb_idx")        );

  int nb_nodes = nodes.extents()[0];

  ArrayView<double,2> latlon ( nodes.field("coordinates") );

  // Identify my master and slave nodes on own partition
  // master nodes are at x=0,  slave nodes are at x=2pi
  std::map<LatLon,int> master_lookup;
  std::map<LatLon,int>  slave_lookup;
  std::vector<int> master_nodes; master_nodes.reserve( 3*nb_nodes );
  std::vector<int>  slave_nodes;  slave_nodes.reserve( 3*nb_nodes );

  int west = 0;
  int east = microdeg( 2.*M_PI );
  for( int jnode=0; jnode<nodes.extents()[0]; ++jnode)
  {
    int x = microdeg(latlon(jnode,XX));
    if( x == west && (part(jnode)==mypart && jnode==loc_idx(jnode))  )
    {
      master_lookup[ LatLon(latlon[jnode]) ] = jnode;
      master_nodes.push_back( x );
      master_nodes.push_back( microdeg(latlon(jnode,YY)) );
      master_nodes.push_back( jnode );
    }
    if( x >= east )
    {
      slave_lookup[ LatLon(latlon[jnode]) ] = jnode;
      slave_nodes.push_back( x );
      slave_nodes.push_back( microdeg(latlon(jnode,YY)) );
      slave_nodes.push_back( jnode );
    }
  }

  std::cout << "found " << master_nodes.size()/3 << " master nodes " << std::endl;
  std::cout << "found " <<  slave_nodes.size()/3 << "  slave nodes " << std::endl;

                      //  std::vector< std::vector<int> > found_slave(MPL::size());
                      //  // Find slaves on other tasks to send to me
                      //  {
                      //    int sendcnt = master_nodes.size();
                      //    std::vector< int > recvcounts( MPL::size() );
                      //    MPL_CHECK_RESULT( MPI_Allgather(&sendcnt,           1, MPI_INT,
                      //                                     recvcounts.data(), 1, MPI_INT, MPI_COMM_WORLD ) );

                      //    std::vector<int> recvdispls( MPL::size() );
                      //    recvdispls[0] = 0;
                      //    int recvcnt = recvcounts[0];
                      //    for( int jproc=1; jproc<MPL::size(); ++jproc )
                      //    {
                      //      recvdispls[jproc] = recvdispls[jproc-1] + recvcounts[jproc-1];
                      //      recvcnt += recvcounts[jproc];
                      //    }
                      //    std::vector<int> recvbuf(recvcnt);

                      //    MPL_CHECK_RESULT( MPI_Allgatherv(
                      //                      master_nodes.data(), master_nodes.size(), MPI_INT,
                      //                      recvbuf.data(), recvcounts.data(), recvdispls.data(),
                      //                      MPI_INT, MPI_COMM_WORLD) );


                      //    for( int jproc=0; jproc<MPL::size(); ++jproc )
                      //    {
                      //      found_slave.reserve(slave_nodes.size());
                      //      ArrayView<int,2> recv_master_latlon(recvbuf.data()+recvdispls[jproc], Extents(recvcounts[jproc]/2,2).data() );
                      //      for( int jnode=0; jnode<recv_master_latlon.extents()[0]; ++jnode )
                      //      {
                      //        LatLon master( recv_master_latlon(jnode,XX)+east, recv_master_latlon(jnode,YY) );
                      //        if( slave_lookup.count( master ) )
                      //        {
                      //          int loc_idx = slave_lookup[ master ];
                      //          found_slave[jproc].push_back( loc_idx );
                      //        }
                      //      }
                      //    }
                      //  }

  std::vector< std::vector<int> > found_master(MPL::size());
  std::vector< std::vector<int> > send_slave_idx(MPL::size());
  // Find masters on other tasks to send to me
  {
    int sendcnt = slave_nodes.size();
    std::vector< int > recvcounts( MPL::size() );
    MPL_CHECK_RESULT( MPI_Allgather(&sendcnt,           1, MPI_INT,
                                     recvcounts.data(), 1, MPI_INT, MPI_COMM_WORLD ) );

    std::vector<int> recvdispls( MPL::size() );
    recvdispls[0] = 0;
    int recvcnt = recvcounts[0];
    for( int jproc=1; jproc<MPL::size(); ++jproc )
    {
      recvdispls[jproc] = recvdispls[jproc-1] + recvcounts[jproc-1];
      recvcnt += recvcounts[jproc];
    }
    std::vector<int> recvbuf(recvcnt);

    MPL_CHECK_RESULT( MPI_Allgatherv(
                      slave_nodes.data(), slave_nodes.size(), MPI_INT,
                      recvbuf.data(), recvcounts.data(), recvdispls.data(),
                      MPI_INT, MPI_COMM_WORLD) );


    for( int jproc=0; jproc<MPL::size(); ++jproc )
    {
      found_master.reserve(master_nodes.size());
      send_slave_idx.reserve(master_nodes.size());
      ArrayView<int,2> recv_slave(recvbuf.data()+recvdispls[jproc], Extents(recvcounts[jproc]/3,3).data() );
      for( int jnode=0; jnode<recv_slave.extents()[0]; ++jnode )
      {
        LatLon slave( recv_slave(jnode,XX)-east, recv_slave(jnode,YY) );
        if( master_lookup.count( slave ) )
        {
          int master_idx = master_lookup[ slave ];
          int slave_idx  = recv_slave(jnode,2);
          found_master[jproc].push_back( master_idx );
          send_slave_idx[jproc].push_back( slave_idx );
          //std::cout << master_idx << " master of " << slave_idx << std::endl;
        }
      }
    }
  }

  // Fill in data to communicate
  std::vector< std::vector<int> > recv_slave_idx( MPL::size() );
  //std::vector< std::vector<int> > send_master_glb_idx( MPL::size() );
  //std::vector< std::vector<int> > recv_master_glb_idx( MPL::size() );
  std::vector< std::vector<int> > send_master_part( MPL::size() );
  std::vector< std::vector<int> > recv_master_part( MPL::size() );
  std::vector< std::vector<int> > send_master_loc( MPL::size() );
  std::vector< std::vector<int> > recv_master_loc( MPL::size() );

                      //  std::vector< std::vector<int> > send_slave_glb_idx( MPL::size() );
                      //  std::vector< std::vector<int> > recv_slave_glb_idx( MPL::size() );
                      //  std::vector< std::vector<int> > send_slave_part( MPL::size() );
                      //  std::vector< std::vector<int> > recv_slave_part( MPL::size() );
                      //  std::vector< std::vector<int> > send_slave_loc( MPL::size() );
                      //  std::vector< std::vector<int> > recv_slave_loc( MPL::size() );

  {
    for( int jproc=0; jproc<MPL::size(); ++jproc )
    {
      int nb_found_master = found_master[jproc].size();
      //send_master_glb_idx[jproc].resize(nb_found_master);
      send_master_part   [jproc].resize(nb_found_master);
      send_master_loc    [jproc].resize(nb_found_master);
      for( int jnode=0; jnode<nb_found_master; ++jnode )
      {
        int loc_idx = found_master[jproc][jnode];
        //send_master_glb_idx[jproc][jnode] = glb_idx( loc_idx );
        send_master_part   [jproc][jnode] = part   ( loc_idx );
        send_master_loc    [jproc][jnode] = loc_idx;
      }

                      //      int nb_found_slaves = found_slave[jproc].size();
                      //      send_slave_glb_idx[jproc].resize(nb_found_slaves);
                      //      send_slave_part   [jproc].resize(nb_found_slaves);
                      //      send_slave_loc    [jproc].resize(nb_found_slaves);
                      //      for( int jnode=0; jnode<nb_found_slaves; ++jnode )
                      //      {
                      //        int loc_idx = found_slave[jproc][jnode];
                      //        send_slave_glb_idx[jproc][jnode] = glb_idx( loc_idx );
                      //        send_slave_part   [jproc][jnode] = part   ( loc_idx );
                      //        send_slave_loc    [jproc][jnode] = loc_idx;
                      //      }
    }
  }

  // Communicate
  MPL::Alltoall( send_slave_idx,      recv_slave_idx      );
  //MPL::Alltoall( send_master_glb_idx, recv_master_glb_idx );
  MPL::Alltoall( send_master_part,    recv_master_part    );
  MPL::Alltoall( send_master_loc,     recv_master_loc     );
                    //  MPL::Alltoall( send_slave_glb_idx,  recv_slave_glb_idx );
                    //  MPL::Alltoall( send_slave_part,     recv_slave_part    );
                    //  MPL::Alltoall( send_slave_loc,      recv_slave_loc     );

  // Fill in periodic
  int nb_recv_master = 0;
  for( int jproc=0; jproc<MPL::size(); ++jproc )
  {
    int nb_recv = recv_slave_idx[jproc].size();
    for( int jnode=0; jnode<nb_recv; ++jnode )
    {
      int slave_idx        = recv_slave_idx     [jproc][jnode];
      //glb_idx( slave_idx ) = recv_master_glb_idx[jproc][jnode];
      part   ( slave_idx ) = recv_master_part   [jproc][jnode];
      loc_idx( slave_idx ) = recv_master_loc    [jproc][jnode];

      //std::cout << part( slave_idx ) <<"["<<loc_idx( slave_idx ) << "] master of " << mypart << "["<<slave_idx<<"]" << std::endl;

    }
  }

}

// ------------------------------------------------------------------

} // namespace actions
} // namespace atlas

