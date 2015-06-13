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

#include "atlas/runtime/ErrorHandling.h"
#include "atlas/mpi/mpi.h"
#include "atlas/Mesh.h"
#include "atlas/FunctionSpace.h"
#include "atlas/actions/BuildPeriodicBoundaries.h"
#include "atlas/Parameters.h"
#include "atlas/util/Array.h"
#include "atlas/util/ArrayView.h"
#include "atlas/util/IndexView.h"
#include "atlas/util/Bitflags.h"
#include "atlas/util/LonLatMicroDeg.h"
#include "atlas/util/PeriodicTransform.h"

using atlas::util::Topology;
using atlas::util::LonLatMicroDeg;
using atlas::util::PeriodicTransform;

namespace atlas {
namespace actions {

typedef gidx_t uid_t;

void build_periodic_boundaries( Mesh& mesh )
{
  int mypart = eckit::mpi::rank();

  FunctionSpace& nodes = mesh.function_space("nodes");

  ArrayView<int,1> flags( nodes.field("flags") );
  IndexView<int,1> ridx ( nodes.field("remote_idx") );
  ArrayView<int,1> part ( nodes.field("partition") );

  int nb_nodes = nodes.shape(0);

  ArrayView<double,2> lonlat ( nodes.field("lonlat") );

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
        LonLatMicroDeg ll(lonlat[jnode]);
        master_lookup[ ll.unique() ] = jnode;
        master_nodes.push_back( ll.lon() );
        master_nodes.push_back( ll.lat() );
        master_nodes.push_back( jnode );
      }
    }
    else if( Topology::check(flags(jnode),Topology::BC|Topology::EAST) )
    {
      Topology::set(flags(jnode),Topology::PERIODIC);
      Topology::set(flags(jnode),Topology::GHOST);
      LonLatMicroDeg ll(lonlat[jnode]);
      slave_lookup[ ll.unique() ] = jnode;
      slave_nodes.push_back( ll.lon() );
      slave_nodes.push_back( ll.lat() );
      slave_nodes.push_back( jnode );
      ridx( jnode ) = -1;
    }
  }

  std::vector< std::vector<int> > found_master(eckit::mpi::size());
  std::vector< std::vector<int> > send_slave_idx(eckit::mpi::size());
  // Find masters on other tasks to send to me
  {
    int sendcnt = slave_nodes.size();
    std::vector< int > recvcounts( eckit::mpi::size() );
    ECKIT_MPI_CHECK_RESULT( MPI_Allgather(&sendcnt,           1, MPI_INT,
                                     recvcounts.data(), 1, MPI_INT, eckit::mpi::comm() ) );

    std::vector<int> recvdispls( eckit::mpi::size() );
    recvdispls[0] = 0;
    int recvcnt = recvcounts[0];
    for( int jproc=1; jproc<eckit::mpi::size(); ++jproc )
    {
      recvdispls[jproc] = recvdispls[jproc-1] + recvcounts[jproc-1];
      recvcnt += recvcounts[jproc];
    }
    std::vector<int> recvbuf(recvcnt);

    ECKIT_MPI_CHECK_RESULT( MPI_Allgatherv(
                      slave_nodes.data(), slave_nodes.size(), eckit::mpi::datatype<int>(),
                      recvbuf.data(), recvcounts.data(), recvdispls.data(),
                      eckit::mpi::datatype<int>(), eckit::mpi::comm()) );


    PeriodicTransform transform;
    for( int jproc=0; jproc<eckit::mpi::size(); ++jproc )
    {
      found_master.reserve(master_nodes.size());
      send_slave_idx.reserve(master_nodes.size());
      ArrayView<int,2> recv_slave(recvbuf.data()+recvdispls[jproc], make_shape(recvcounts[jproc]/3,3).data() );
      for( int jnode=0; jnode<recv_slave.shape(0); ++jnode )
      {
        LonLatMicroDeg slave( recv_slave(jnode,LON), recv_slave(jnode,LAT) );
        transform(slave,-1);
        uid_t slave_uid = slave.unique();
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
  std::vector< std::vector<int> > recv_slave_idx( eckit::mpi::size() );
  std::vector< std::vector<int> > send_master_part( eckit::mpi::size() );
  std::vector< std::vector<int> > recv_master_part( eckit::mpi::size() );
  std::vector< std::vector<int> > send_master_ridx( eckit::mpi::size() );
  std::vector< std::vector<int> > recv_master_ridx( eckit::mpi::size() );

                      //  std::vector< std::vector<int> > send_slave_part( eckit::mpi::size() );
                      //  std::vector< std::vector<int> > recv_slave_part( eckit::mpi::size() );
                      //  std::vector< std::vector<int> > send_slave_ridx( eckit::mpi::size() );
                      //  std::vector< std::vector<int> > recv_slave_ridx( eckit::mpi::size() );

  {
    for( int jproc=0; jproc<eckit::mpi::size(); ++jproc )
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
  eckit::mpi::all_to_all( send_slave_idx,      recv_slave_idx      );
  eckit::mpi::all_to_all( send_master_part,    recv_master_part    );
  eckit::mpi::all_to_all( send_master_ridx,    recv_master_ridx     );
                    //  eckit::mpi::all_to_all( send_slave_part,     recv_slave_part    );
                    //  eckit::mpi::all_to_all( send_slave_loc,      recv_slave_ridx    );

  // Fill in periodic
  int nb_recv_master = 0;
  for( int jproc=0; jproc<eckit::mpi::size(); ++jproc )
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
  ATLAS_ERROR_HANDLING( build_periodic_boundaries(*mesh) );
}
// ------------------------------------------------------------------



} // namespace actions
} // namespace atlas

