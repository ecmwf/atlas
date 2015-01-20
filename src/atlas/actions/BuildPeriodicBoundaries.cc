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
#include <cmath>

#include "atlas/Mesh.h"
#include "atlas/FunctionSpace.h"
#include "atlas/actions/BuildPeriodicBoundaries.h"
#include "atlas/Parameters.h"
#include "atlas/Util.h"
#include "atlas/util/Array.h"
#include "atlas/util/ArrayView.h"
#include "atlas/util/IndexView.h"

namespace atlas {
namespace actions {

typedef gidx_t uid_t;

void build_periodic_boundaries( Mesh& mesh )
{
  int mypart = mpi::rank();

  FunctionSpace& nodes = mesh.function_space("nodes");

  ArrayView<int,1> flags( nodes.field("flags") );
  IndexView<int,1> ridx ( nodes.field("remote_idx") );
  ArrayView<int,1> part ( nodes.field("partition") );

  int nb_nodes = nodes.shape(0);

  ArrayView<double,2> latlon ( nodes.field("coordinates") );

  // Identify my master and slave nodes on own partition
  // master nodes are at x=0,  slave nodes are at x=2pi
  std::map<uid_t,int> master_lookup;
  std::map<uid_t,int>  slave_lookup;
  std::vector<int> master_nodes; master_nodes.reserve( 3*nb_nodes );
  std::vector<int> slave_nodes;  slave_nodes.reserve( 3*nb_nodes );

  for( int jnode=0; jnode<nodes.shape(0); ++jnode)
  {
    if( Topology::check_all(flags(jnode),Topology::BC|Topology::WEST) )
    {
      if( part(jnode) == mypart )
      {
        Topology::set(flags(jnode),Topology::PERIODIC);
        LatLonPoint ll(latlon[jnode]);
        master_lookup[ ll.uid() ] = jnode;
        master_nodes.push_back( ll.x );
        master_nodes.push_back( ll.y );
        master_nodes.push_back( jnode );
      }
    }
    else if( Topology::check(flags(jnode),Topology::BC|Topology::EAST) )
    {
      Topology::set(flags(jnode),Topology::PERIODIC);
      Topology::set(flags(jnode),Topology::GHOST);
      LatLonPoint ll(latlon[jnode]);
      slave_lookup[ ll.uid() ] = jnode;
      slave_nodes.push_back( ll.x );
      slave_nodes.push_back( ll.y );
      slave_nodes.push_back( jnode );
      ridx( jnode ) = -1;
    }
  }

  std::vector< std::vector<int> > found_master(mpi::size());
  std::vector< std::vector<int> > send_slave_idx(mpi::size());
  // Find masters on other tasks to send to me
  {
    int sendcnt = slave_nodes.size();
    std::vector< int > recvcounts( mpi::size() );
    MPL_CHECK_RESULT( MPI_Allgather(&sendcnt,           1, MPI_INT,
                                     recvcounts.data(), 1, MPI_INT, MPI_COMM_WORLD ) );

    std::vector<int> recvdispls( mpi::size() );
    recvdispls[0] = 0;
    int recvcnt = recvcounts[0];
    for( int jproc=1; jproc<mpi::size(); ++jproc )
    {
      recvdispls[jproc] = recvdispls[jproc-1] + recvcounts[jproc-1];
      recvcnt += recvcounts[jproc];
    }
    std::vector<int> recvbuf(recvcnt);

    MPL_CHECK_RESULT( MPI_Allgatherv(
                      slave_nodes.data(), slave_nodes.size(), mpi::TYPE<int>(),
                      recvbuf.data(), recvcounts.data(), recvdispls.data(),
                      mpi::TYPE<int>(), MPI_COMM_WORLD) );


    PeriodicTransform transform;
    for( int jproc=0; jproc<mpi::size(); ++jproc )
    {
      found_master.reserve(master_nodes.size());
      send_slave_idx.reserve(master_nodes.size());
      ArrayView<int,2> recv_slave(recvbuf.data()+recvdispls[jproc], make_shape(recvcounts[jproc]/3,3).data() );
      for( int jnode=0; jnode<recv_slave.shape(0); ++jnode )
      {
        LatLonPoint slave( recv_slave(jnode,XX), recv_slave(jnode,YY) );
        transform(slave,-1);
        uid_t slave_uid = slave.uid();
        if( master_lookup.count( slave_uid ) )
        {
          int master_idx = master_lookup[ slave_uid ];
          int slave_idx  = recv_slave(jnode,2);
          found_master[jproc].push_back( master_idx );
          send_slave_idx[jproc].push_back( slave_idx );
        }
      }
    }
  }

  // Fill in data to communicate
  std::vector< std::vector<int> > recv_slave_idx( mpi::size() );
  std::vector< std::vector<int> > send_master_part( mpi::size() );
  std::vector< std::vector<int> > recv_master_part( mpi::size() );
  std::vector< std::vector<int> > send_master_ridx( mpi::size() );
  std::vector< std::vector<int> > recv_master_ridx( mpi::size() );

                      //  std::vector< std::vector<int> > send_slave_part( mpi::size() );
                      //  std::vector< std::vector<int> > recv_slave_part( mpi::size() );
                      //  std::vector< std::vector<int> > send_slave_ridx( mpi::size() );
                      //  std::vector< std::vector<int> > recv_slave_ridx( mpi::size() );

  {
    for( int jproc=0; jproc<mpi::size(); ++jproc )
    {
      int nb_found_master = found_master[jproc].size();
      send_master_part   [jproc].resize(nb_found_master);
      send_master_ridx    [jproc].resize(nb_found_master);
      for( int jnode=0; jnode<nb_found_master; ++jnode )
      {
        int loc_idx = found_master[jproc][jnode];
        send_master_part[jproc][jnode] = part   ( loc_idx );
        send_master_ridx[jproc][jnode] = loc_idx;
      }

                      //      int nb_found_slaves = found_slave[jproc].size();
                      //      send_slave_glb_idx[jproc].resize(nb_found_slaves);
                      //      send_slave_part   [jproc].resize(nb_found_slaves);
                      //      send_slave_ridx   [jproc].resize(nb_found_slaves);
                      //      for( int jnode=0; jnode<nb_found_slaves; ++jnode )
                      //      {
                      //        int loc_idx = found_slave[jproc][jnode];
                      //        send_slave_glb_idx[jproc][jnode] = glb_idx( loc_idx );
                      //        send_slave_part   [jproc][jnode] = part   ( loc_idx );
                      //        send_slave_ridx   [jproc][jnode] = loc_idx;
                      //      }
    }
  }

  // Communicate
  mpi::Alltoall( send_slave_idx,      recv_slave_idx      );
  mpi::Alltoall( send_master_part,    recv_master_part    );
  mpi::Alltoall( send_master_ridx,    recv_master_ridx     );
                    //  mpi::Alltoall( send_slave_part,     recv_slave_part    );
                    //  mpi::Alltoall( send_slave_loc,      recv_slave_ridx    );

  // Fill in periodic
  int nb_recv_master = 0;
  for( int jproc=0; jproc<mpi::size(); ++jproc )
  {
    int nb_recv = recv_slave_idx[jproc].size();
    for( int jnode=0; jnode<nb_recv; ++jnode )
    {
      int slave_idx     = recv_slave_idx  [jproc][jnode];
      part( slave_idx ) = recv_master_part[jproc][jnode];
      ridx( slave_idx ) = recv_master_ridx[jproc][jnode];
    }
  }

}

// ------------------------------------------------------------------
// C wrapper interfaces to C++ routines

void atlas__build_periodic_boundaries ( Mesh* mesh) {
  build_periodic_boundaries(*mesh);
}
// ------------------------------------------------------------------



} // namespace actions
} // namespace atlas

