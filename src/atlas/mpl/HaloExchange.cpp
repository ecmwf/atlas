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

#include "atlas/atlas_defines.h"
#include "atlas/mpl/HaloExchange.hpp"

#ifdef HAVE_FORTRAN_NUMBERING
#define FORTRAN 1+
#else
#define FORTRAN
#endif

namespace atlas {

namespace {
struct IsGhostPoint
{
  IsGhostPoint( const int part[], const int ridx[], const int N )
  {
    part_   = part;
    ridx_   = ridx;
    mypart_ = MPL::rank();
  }

  bool operator()(int idx)
  {
    if( part_[idx] != mypart_     ) return true;
    if( ridx_[idx] != FORTRAN idx ) return true;
    return false;
  }
  int mypart_;
  const int* part_;
  const int* ridx_;
};
}

HaloExchange::HaloExchange() :
  is_setup_(false)
{
  myproc = MPL::rank();
  nproc  = MPL::size();
}

void HaloExchange::setup( const int part[],
                          const int remote_idx[],
                          const int size )
{
  parsize_ = size;
  sendcounts_.resize(nproc,0);
  recvcounts_.resize(nproc,0);
  senddispls_.resize(nproc,0);
  recvdispls_.resize(nproc,0);

  /*
    Find the amount of nodes this proc has to receive from each other proc
  */

  IsGhostPoint is_ghost(part,remote_idx,size);

  for (int jj=0; jj<size; ++jj)
  {
    if ( is_ghost(jj) )
      ++recvcounts_[part[jj]];
  }
  recvcnt_ = std::accumulate(recvcounts_.begin(),recvcounts_.end(),0);
//  std::cout << myproc << ":  recvcnt = " << recvcnt_ << std::endl;


  /*
    Find the amount of nodes this proc has to send to each other proc
  */

  MPL_CHECK_RESULT( MPI_Alltoall( recvcounts_.data(), 1, MPI_INT, sendcounts_.data(), 1, MPI_INT, MPI_COMM_WORLD ) );
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
    Fill vector "send_requests" with remote index of nodes needed, but are on other procs
    We can also fill in the vector "recvmap_" which holds local indices of requested nodes
  */

  std::vector<int> send_requests(recvcnt_);

  recvmap_.resize(recvcnt_);
  std::vector<int> cnt(nproc,0);
  for (int jj=0; jj<size; ++jj)
  {
    if ( is_ghost(jj) )
    {
      const int req_idx = recvdispls_[part[jj]] + cnt[part[jj]];
      send_requests[req_idx] = remote_idx[jj];
      recvmap_[req_idx] = jj;
      ++cnt[part[jj]];
    }
  }

  /*
    Fill vector "recv_requests" with what is needed by other procs
  */

  std::vector<int> recv_requests(sendcnt_);

  MPL_CHECK_RESULT( MPI_Alltoallv(
                      send_requests.data(), recvcounts_.data(), recvdispls_.data(), MPI_INT,
                      recv_requests.data(), sendcounts_.data(), senddispls_.data(), MPI_INT,
                      MPI_COMM_WORLD ) );

  /*
    What needs to be sent to other procs is asked by remote_idx, which is local here
  */
  sendmap_.resize(sendcnt_);
  for( int jj=0; jj<sendcnt_; ++jj )
    sendmap_[jj] = recv_requests[jj];

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



template <>
void HaloExchange::execute<int,1>( ArrayView<int,1>& field ) const
{
  int one[] = {1};
  execute( field.data(), one, one, 1 );
}

template <>
void HaloExchange::execute<float,1>( ArrayView<float,1>& field ) const
{
  int one[] = {1};
  execute( field.data(), one, one, 1 );
}

template <>
void HaloExchange::execute<double,1>( ArrayView<double,1>& field ) const
{
  int one[] = {1};
  execute( field.data(), one, one, 1 );
}

/////////////////////


HaloExchange* atlas__HaloExchange__new () { 
  return new HaloExchange(); 
}

void atlas__HaloExchange__delete (HaloExchange* This) {
  delete This;
}

void atlas__HaloExchange__setup (HaloExchange* This, int part[], int remote_idx[], int size)
{
  This->setup(part,remote_idx,size);
}

void atlas__HaloExchange__execute_strided_int (HaloExchange* This, int field[], int var_strides[], int var_extents[], int var_rank) {
  This->execute(field,var_strides,var_extents,var_rank);
}

void atlas__HaloExchange__execute_strided_double (HaloExchange* This, float field[], int var_strides[], int var_extents[], int var_rank) {
  This->execute(field,var_strides,var_extents,var_rank);
}

void atlas__HaloExchange__execute_strided_float (HaloExchange* This, double field[], int var_strides[], int var_extents[], int var_rank) {
  This->execute(field,var_strides,var_extents,var_rank);
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
