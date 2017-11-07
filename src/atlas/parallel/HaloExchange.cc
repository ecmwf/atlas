/*
 * (C) Copyright 1996-2017 ECMWF.
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
#include <sstream>
#include "atlas/parallel/HaloExchange.h"
#include "atlas/parallel/mpi/Statistics.h"
#include "atlas/array/Array.h"

namespace atlas {
namespace parallel {

namespace {
struct IsGhostPoint
{
  IsGhostPoint( const int part[], const int ridx[], const int base, const int N )
  {
    part_   = part;
    ridx_   = ridx;
    base_   = base;
    mypart_ = parallel::mpi::comm().rank();
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
  myproc = parallel::mpi::comm().rank();
  nproc  = parallel::mpi::comm().size();
}

HaloExchange::HaloExchange(const std::string& name) :
  name_(name),
  is_setup_(false)
{
  myproc = parallel::mpi::comm().rank();
  nproc  = parallel::mpi::comm().size();
}

void HaloExchange::setup( const int part[],
                          const int remote_idx[], const int base,
                          const size_t parsize )
{
  ATLAS_TRACE("HaloExchange::setup");

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
  ATLAS_TRACE_MPI( ALLTOALL ) {
    parallel::mpi::comm().allToAll(recvcounts_, sendcounts_);
  }

  sendcnt_ = std::accumulate(sendcounts_.begin(),sendcounts_.end(),0);

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
  ATLAS_TRACE_MPI( ALLTOALL ) {
    parallel::mpi::comm().allToAllv(send_requests.data(), recvcounts_.data(), recvdispls_.data(),
                                    recv_requests.data(), sendcounts_.data(), senddispls_.data());
  }

  /*
    What needs to be sent to other procs is asked by remote_idx, which is local here
  */
  sendmap_.resize(sendcnt_);
  for( int jj=0; jj<sendcnt_; ++jj )
    sendmap_[jj] = recv_requests[jj];

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
    array::ArrayShape shape;
    shape.resize(var_rank);
    shape.assign(var_extents, var_extents+var_rank);

    array::ArrayStrides strides;
    strides.resize(var_rank);
    strides.assign(var_strides, var_strides+var_rank);

    eckit::SharedPtr<array::Array> arr ( array::Array::wrap<int>(field,
        array::ArraySpec{shape, strides}) );

    switch(var_rank) {
        case 1: {This->execute<int,1>(*arr); break;}
        case 2: {This->execute<int,2>(*arr); break;}
        case 3: {This->execute<int,3>(*arr); break;}
        case 4: {This->execute<int,4>(*arr); break;}
        default: throw eckit::AssertionFailed("Rank not supported in halo exchange");
    }
}

void atlas__HaloExchange__execute_strided_long (HaloExchange* This, long field[], int var_strides[], int var_extents[], int var_rank) {
    array::ArrayShape shape;
    shape.resize(var_rank);
    shape.assign(var_extents, var_extents+var_rank);

    array::ArrayStrides strides;
    strides.resize(var_rank);
    strides.assign(var_strides, var_strides+var_rank);

    eckit::SharedPtr<array::Array> arr ( array::Array::wrap<long>(field,
        array::ArraySpec{shape, strides}) );

    switch(var_rank) {
        case 1: {This->execute<long,1>(*arr); break;}
        case 2: {This->execute<long,2>(*arr); break;}
        case 3: {This->execute<long,3>(*arr); break;}
        case 4: {This->execute<long,4>(*arr); break;}
        default: throw eckit::AssertionFailed("Rank not supported in halo exchange");
    }
}

void atlas__HaloExchange__execute_strided_float (HaloExchange* This, float field[], int var_strides[], int var_extents[], int var_rank) {
    array::ArrayShape shape;
    shape.resize(var_rank);
    shape.assign(var_extents, var_extents+var_rank);

    array::ArrayStrides strides;
    strides.resize(var_rank);
    strides.assign(var_strides, var_strides+var_rank);

    eckit::SharedPtr<array::Array> arr ( array::Array::wrap<float>(field,
        array::ArraySpec{shape, strides}) );

    switch(var_rank) {
        case 1: {This->execute<float,1>(*arr); break;}
        case 2: {This->execute<float,2>(*arr); break;}
        case 3: {This->execute<float,3>(*arr); break;}
        case 4: {This->execute<float,4>(*arr); break;}
        default: throw eckit::AssertionFailed("Rank not supported in halo exchange");
    }
}

void atlas__HaloExchange__execute_strided_double (HaloExchange* This, double field[], int var_strides[], int var_extents[], int var_rank) {
    array::ArrayShape shape;
    shape.resize(var_rank);
    shape.assign(var_extents, var_extents+var_rank);

    array::ArrayStrides strides;
    strides.resize(var_rank);
    strides.assign(var_strides, var_strides+var_rank);

    eckit::SharedPtr<array::Array> arr ( array::Array::wrap<double>(field,
        array::ArraySpec{shape, strides}) );

    switch(var_rank) {
        case 1: {This->execute<double,1>(*arr); break;}
        case 2: {This->execute<double,2>(*arr); break;}
        case 3: {This->execute<double,3>(*arr); break;}
        case 4: {This->execute<double,4>(*arr); break;}
        default: throw eckit::AssertionFailed("Rank not supported in halo exchange");
    }
}

/////////////////////

} // namespace parallel
} // namespace atlas
