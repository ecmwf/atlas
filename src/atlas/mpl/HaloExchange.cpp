/*
 * (C) Copyright 1996-2014 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */



#include <stdexcept>
#include <numeric>      // std::accumulate
#include <iostream>
#include <sstream>

#include "atlas/mpl/HaloExchange.hpp"

namespace atlas {

HaloExchange::HaloExchange() :
  is_setup_(false)
{
  myproc = MPL::rank();
  nproc  = MPL::size();
}

void HaloExchange::setup(const int proc[],
                         const int glb_idx[],
                         const int master_glb_idx[],
                         const std::vector<int>& bounds,
                         int par_bound )
{
//  bounds_.resize(bounds.size());
//  for(int i=0; i<bounds.size(); ++i)
//    bounds_[i] = bounds[i];

  bounds_ = bounds;
  par_bound_ = par_bound;


  sendcounts_.resize(nproc,0);
  recvcounts_.resize(nproc,0);
  senddispls_.resize(nproc,0);
  recvdispls_.resize(nproc,0);

  int nb_nodes = bounds_[par_bound_];

  /*
    Create a temporary mapping from global to local indices
    Currently this is quickly implemented using a LONG list...
  */

  int max_glb_idx = -1;
  for (int jj=0; jj<nb_nodes; ++jj)
  {
    max_glb_idx = std::max( max_glb_idx, glb_idx[jj] );
  }

  MPL_CHECK_RESULT( MPI_Allreduce( MPI_IN_PLACE, &max_glb_idx, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD ) );

  std::vector<int> map_glb_to_loc(max_glb_idx+1,-1);
  for (int jj=0; jj<nb_nodes; ++jj)
    map_glb_to_loc[glb_idx[jj]] = jj;

//   std::cout << myproc << ":  map_glb_to_loc = ";
//   for (int i=0; i<max_glb_idx+1; ++i)
//     std::cout << map_glb_to_loc[i] << " ";
//   std::cout << std::endl;

  /*
    Find the amount of nodes this proc has to receive from each other proc
  */

  for (int jj=0; jj<nb_nodes; ++jj)
  {
    if (proc[jj] != myproc || master_glb_idx[jj] != glb_idx[jj])
      ++recvcounts_[proc[jj]];
  }
  recvcnt_ = std::accumulate(recvcounts_.begin(),recvcounts_.end(),0);
//  std::cout << myproc << ":  recvcnt = " << recvcnt_ << std::endl;


  /*
    Find the amount of nodes this proc has to send to each other proc
  */

  MPL_CHECK_RESULT( MPI_Alltoall( &recvcounts_[0], 1, MPI_INT, &sendcounts_[0], 1, MPI_INT, MPI_COMM_WORLD ) );
  sendcnt_ = std::accumulate(sendcounts_.begin(),sendcounts_.end(),0);
  //std::cout << myproc << ":  sendcnt = " << sendcnt_ << std::endl;

  recvdispls_[0]=0;
  senddispls_[0]=0;
  for (int jproc=1; jproc<nproc; ++jproc) // start at 1
  {
    recvdispls_[jproc]=recvcounts_[jproc-1]+recvdispls_[jproc-1];
    senddispls_[jproc]=sendcounts_[jproc-1]+senddispls_[jproc-1];
  }
  /*
    Fill vector "send_requests" with global indices of nodes needed, but are on other procs
    We can also fill in the vector "recvmap_" which holds local indices of requested nodes
  */

  std::vector<int> send_requests(recvcnt_);

  recvmap_.resize(recvcnt_);
  std::vector<int> cnt(nproc,0);
  for (int jj=0; jj<nb_nodes; ++jj)
  {
    if (proc[jj] != myproc || master_glb_idx[jj] != glb_idx[jj])
    {
      const int req_idx = recvdispls_[proc[jj]] + cnt[proc[jj]];
      send_requests[req_idx] = master_glb_idx[jj];
      recvmap_[req_idx] = jj;
      ++cnt[proc[jj]];
    }
//    if (master_glb_idx[jj] != glb_idx[jj])
//    {
//      const int req_idx = recvdispls_[proc[jj]] + cnt[proc[jj]];
//      send_requests[req_idx] = master_glb_idx[jj];
//      recvmap_[req_idx] = jj;
//      ++cnt[proc[jj]];
//    }

  }

  /*
    Fill vector "recv_requests" with global indices that are needed by other procs
  */

  std::vector<int> recv_requests(sendcnt_);

  MPL_CHECK_RESULT( MPI_Alltoallv( &send_requests[0], &recvcounts_[0], &recvdispls_[0], MPI_INT,
                        &recv_requests[0], &sendcounts_[0], &senddispls_[0], MPI_INT,
                        MPI_COMM_WORLD ) );

  /*
    What needs to be sent to other procs can be found by a map from global to local indices
  */
  sendmap_.resize(sendcnt_);
  for( int jj=0; jj<sendcnt_; ++jj )
    sendmap_[jj] = map_glb_to_loc[ recv_requests[jj] ];

  // Packet size
  packet_size_ = 1;
  const int nb_bounds = bounds_.size();
  for (int b=0; b<nb_bounds; ++b)
  {
    if ( b != par_bound_ && bounds_[b] >= 0)
      packet_size_ *= bounds_[b];
  }

//   std::cout << myproc << "  :  sendmap_  = ";
//   for( int i=0; i< sendmap_.size(); ++i)
//     std::cout << sendmap_[i] << " ";
//   std::cout << std::endl;

//   std::cout << myproc << "  :  recvmap_  = ";
//   for( int i=0; i< recvmap_.size(); ++i)
//     std::cout << recvmap_[i] << " ";
//   std::cout << std::endl;

  is_setup_ = true;
}



template<>
inline void HaloExchange::create_mappings_impl<2,0>( 
    std::vector<int>& send_map, 
    std::vector<int>& recv_map,
    int nb_vars) const
{
  const int nb_nodes = bounds_[0];
  int send_idx(0);
  for (int jnode=0; jnode<sendcnt_; ++jnode)
  {
    const int inode = sendmap_[jnode];
    for (int jvar=0; jvar<nb_vars; ++jvar)
    {
      send_map[send_idx++] = index( inode,jvar,   nb_nodes,nb_vars);
    }
  }
  int recv_idx(0);
  for (int jnode=0; jnode<recvcnt_; ++jnode)
  {
    const int inode = recvmap_[jnode];
    for (int jvar=0; jvar<nb_vars; ++jvar)
    {
      recv_map[recv_idx++] = index( inode,jvar,   nb_nodes,nb_vars);;
    }
  }
}

template<>
inline void HaloExchange::create_mappings_impl<2,1>( 
    std::vector<int>& send_map, 
    std::vector<int>& recv_map,
    int nb_vars) const
{
  const int nb_nodes = bounds_[1];
  int send_idx(0);
  for (int jnode=0; jnode<sendcnt_; ++jnode)
  {
    const int inode = sendmap_[jnode];
    for (int jvar=0; jvar<nb_vars; ++jvar)
    {
      send_map[send_idx++] = index( jvar, inode,  nb_vars, nb_nodes);
    }
  }
  int recv_idx(0);
  for (int jnode=0; jnode<recvcnt_; ++jnode)
  {
    const int inode = recvmap_[jnode];
    for (int jvar=0; jvar<nb_vars; ++jvar)
    {
      recv_map[recv_idx++] = index( jvar, inode,  nb_vars, nb_nodes);
    }
  }
}

/// create_mappings_impl<3,1>
template<>
inline void HaloExchange::create_mappings_impl<3,1>( 
    std::vector<int>& send_map, 
    std::vector<int>& recv_map,
    int nb_vars) const
{
  const int nb_levs = bounds_[0];
  const int nb_nodes = bounds_[1];
  int send_idx(0);
  int recv_idx(0);
  for (int n=0; n<sendcnt_; ++n)
  {
    const int jnode = sendmap_[n];
    for (int var=0; var<nb_vars; ++var)
    {
      const int varidx = var*nb_nodes;
      const int nodevaridx = (jnode + varidx)*nb_levs;
      for (int l=0; l<nb_levs; ++l)
      {
        const int levnodevaridx = l + nodevaridx;
        send_map[send_idx++] = levnodevaridx;
      }
    }
  }

  for (int n=0; n<recvcnt_; ++n)
  {
    const int jnode = recvmap_[n];
    for (int var=0; var<nb_vars; ++var)
    {
      const int varidx = var*nb_nodes;
      const int nodevaridx = (jnode + varidx)*nb_levs;
      for (int l=0; l<nb_levs; ++l)
      {
        const int levnodevaridx = l + nodevaridx;
        recv_map[recv_idx++] = levnodevaridx;
      }
    }
  }
}

void HaloExchange::create_mappings( 
    std::vector<int>& send_map, 
    std::vector<int>& recv_map,
    int nb_vars ) const
{
  const int nb_bounds = bounds_.size();
  switch (nb_bounds)
  {
    case 2:
    {
      switch (par_bound_)
      {
        case 0:
          create_mappings_impl<2,0>(send_map,recv_map,nb_vars);
          break;
        case 1:
          create_mappings_impl<2,1>(send_map,recv_map,nb_vars);
          break;
        default:
          std::stringstream errmsg;
          errmsg << "create_mappings<"<<nb_bounds<<","<<par_bound_<<"> not implemented";
          throw std::runtime_error(errmsg.str());
      }
      break;
    }
    case 3:
    {
      switch (par_bound_)
      {
        case 1:
          create_mappings_impl<3,1>(send_map,recv_map,nb_vars);
          break;
        default:
          std::stringstream errmsg;
          errmsg << "create_mappings<"<<nb_bounds<<","<<par_bound_<<"> not implemented";
          throw std::runtime_error(errmsg.str());
      }
      break;
    }
    default:
      std::stringstream errmsg;
      errmsg << "create_mappings<"<<nb_bounds<<","<<par_bound_<<"> not implemented";
      throw std::runtime_error(errmsg.str());
  }
}

/////////////////////


HaloExchange* atlas__HaloExchange__new () { 
  return new HaloExchange(); 
}

void atlas__HaloExchange__delete (HaloExchange* This) {
  delete This;
}

void atlas__HaloExchange__setup (HaloExchange* This, int proc[], int glb_idx[], int master_glb_idx[], int bounds[], int nb_bounds, int par_bound)
{
  std::vector<int> bounds_vec(bounds,bounds+nb_bounds);
  This->setup(proc,glb_idx,master_glb_idx,bounds_vec,par_bound);
}

void atlas__HaloExchange__execute_int (HaloExchange* This, int field[], int nb_vars ) { 
  This->execute(field,nb_vars);
}

void atlas__HaloExchange__execute_float (HaloExchange* This, float field[], int nb_vars ) { 
  This->execute(field,nb_vars);
}

void atlas__HaloExchange__execute_double (HaloExchange* This, double field[], int nb_vars ) { 
  This->execute(field,nb_vars);
}

/////////////////////

}
