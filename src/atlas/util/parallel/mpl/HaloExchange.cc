/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */


/// @author Willem Deconinck
/// @date   Nov 2013

#include <stdexcept>
#include <numeric>
#include <iostream>
#include <sstream>
#include "atlas/util/parallel/mpl/HaloExchange.h"

namespace atlas {
namespace mpl {
namespace {
struct IsGhostPoint
{
  IsGhostPoint( const int part[], const int ridx[], const int base, const int N )
  {
    part_   = part;
    ridx_   = ridx;
    base_   = base;
    mypart_ = eckit::mpi::rank();
  }

  bool operator()(size_t idx)
  {
    if( part_[idx] != mypart_  ) return true;
    if( size_t(ridx_[idx]) != base_+idx ) return true;
    return false;
  }
  int mypart_;
  const int* part_;
  const int* ridx_;
  int base_;
};
}

HaloExchange::HaloExchange() :
  name_(),
  is_setup_(false)
{
  myproc = eckit::mpi::rank();
  nproc  = eckit::mpi::size();
}

HaloExchange::HaloExchange(const std::string& name) :
  name_(name),
  is_setup_(false)
{
  myproc = eckit::mpi::rank();
  nproc  = eckit::mpi::size();
}

void HaloExchange::setup( const int part[],
                          const int remote_idx[], const int base,
                          const size_t parsize )
{
  parsize_ = parsize;
  sendcounts_.resize(nproc); sendcounts_.assign(nproc,0);
  recvcounts_.resize(nproc); recvcounts_.assign(nproc,0);
  senddispls_.resize(nproc); senddispls_.assign(nproc,0);
  recvdispls_.resize(nproc); recvdispls_.assign(nproc,0);

  /*
    Find the amount of nodes this proc has to receive from each other proc
  */

  IsGhostPoint is_ghost(part,remote_idx,base,parsize_);

  for (int jj = 0; jj < parsize_; ++jj)
  {
    if ( is_ghost(jj) )
      ++recvcounts_[part[jj]];
  }
  recvcnt_ = std::accumulate(recvcounts_.begin(),recvcounts_.end(),0);


  /*
    Find the amount of nodes this proc has to send to each other proc
  */

  ECKIT_MPI_CHECK_RESULT( MPI_Alltoall( recvcounts_.data(), 1, MPI_INT, sendcounts_.data(), 1, MPI_INT, eckit::mpi::comm() ) );
  sendcnt_ = std::accumulate(sendcounts_.begin(),sendcounts_.end(),0);
//  std::cout << myproc << ":  sendcnt = " << sendcnt_ << std::endl;
//  std::cout << myproc << ":  recvcnt = " << recvcnt_ << std::endl;

  recvdispls_[0]=0;
  senddispls_[0]=0;
  for (int jproc = 1; jproc < nproc; ++jproc) // start at 1
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
  for (int jj = 0; jj < parsize_; ++jj)
  {
    if ( is_ghost(jj) )
    {
      const int req_idx = recvdispls_[part[jj]] + cnt[part[jj]];
      send_requests[req_idx] = remote_idx[jj]-base;
      recvmap_[req_idx] = jj;
      ++cnt[part[jj]];
    }
  }

  /*
    Fill vector "recv_requests" with what is needed by other procs
  */

  std::vector<int> recv_requests(sendcnt_);

  ECKIT_MPI_CHECK_RESULT( MPI_Alltoallv(
                      send_requests.data(), recvcounts_.data(), recvdispls_.data(), MPI_INT,
                      recv_requests.data(), sendcounts_.data(), senddispls_.data(), MPI_INT,
                      eckit::mpi::comm() ) );

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


/////////////////////


HaloExchange* atlas__HaloExchange__new () {
  return new HaloExchange();
}

void atlas__HaloExchange__delete (HaloExchange* This) {
  delete This;
}

void atlas__HaloExchange__setup (HaloExchange* This, int part[], int remote_idx[], int base, int size)
{
  This->setup(part,remote_idx,base,size);
}

void atlas__HaloExchange__execute_strided_int (HaloExchange* This, int field[], int var_strides[], int var_extents[], int var_rank) {
  std::vector<size_t> vstrides(var_rank);
  std::vector<size_t> vextents(var_rank);
  for(int n = 0; n < var_rank; ++n) {
    vstrides[n] = var_strides[n];
    vextents[n] = var_extents[n];
  }
  This->execute(field,vstrides.data(),vextents.data(),var_rank);
}

void atlas__HaloExchange__execute_strided_long (HaloExchange* This, long field[], int var_strides[], int var_extents[], int var_rank) {
  std::vector<size_t> vstrides(var_rank);
  std::vector<size_t> vextents(var_rank);
  for(int n = 0; n < var_rank; ++n) {
    vstrides[n] = var_strides[n];
    vextents[n] = var_extents[n];
  }
  This->execute(field,vstrides.data(),vextents.data(),var_rank);
}

void atlas__HaloExchange__execute_strided_float (HaloExchange* This, float field[], int var_strides[], int var_extents[], int var_rank) {
  std::vector<size_t> vstrides(var_rank);
  std::vector<size_t> vextents(var_rank);
  for(int n = 0; n < var_rank; ++n) {
    vstrides[n] = var_strides[n];
    vextents[n] = var_extents[n];
  }
  This->execute(field,vstrides.data(),vextents.data(),var_rank);
}

void atlas__HaloExchange__execute_strided_double (HaloExchange* This, double field[], int var_strides[], int var_extents[], int var_rank) {
  std::vector<size_t> vstrides(var_rank);
  std::vector<size_t> vextents(var_rank);
  for(int n = 0; n < var_rank; ++n) {
    vstrides[n] = var_strides[n];
    vextents[n] = var_extents[n];
  }
  This->execute(field,vstrides.data(),vextents.data(),var_rank);
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

} // namespace mpl
} // namespace atlas
