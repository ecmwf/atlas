/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <stdexcept>
#include <numeric>
#include <iostream>
#include <sstream>
#include "atlas/internals/Checksum.h"
#include "atlas/array/Array.h"
#include "atlas/array/ArrayView.h"
#include "atlas/runtime/Log.h"
#include "atlas/parallel/GatherScatter.h"

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

  bool operator()(int idx)
  {
    if( part_[idx] != mypart_  ) return true;
    if( ridx_[idx] != base_+idx ) return true;
    return false;
  }
  int mypart_;
  const int* part_;
  const int* ridx_;
  int base_;
};

struct Node
{
  gidx_t p,i;
  gidx_t g;

  Node() {}
  Node(gidx_t gid, int part, int idx)
  {
    g = gid;
    p = part;
    i = idx;
  }

  bool operator < (const Node& other) const
  {
    return ( g < other.g );
  }

  bool operator == (const Node& other) const
  {
    return ( g == other.g );
  }
};


bool operator < (const gidx_t g, const Node& n)
{
  return ( g < n.g );
}

}

GatherScatter::GatherScatter() :
  name_(),
  is_setup_(false)
{
  myproc = parallel::mpi::comm().rank();
  nproc  = parallel::mpi::comm().size();
}

GatherScatter::GatherScatter(const std::string& name) :
  name_(name),
  is_setup_(false)
{
  myproc = parallel::mpi::comm().rank();
  nproc  = parallel::mpi::comm().size();
}


void GatherScatter::setup( const int part[],
                           const int remote_idx[], const int base,
                           const gidx_t glb_idx[], const int mask[], const size_t parsize )
{
  parsize_ = parsize;

  glbcounts_.resize(nproc); glbcounts_.assign(nproc,0);
  glbdispls_.resize(nproc); glbdispls_.assign(nproc,0);
  const size_t nvar = 3;

  std::vector<gidx_t> sendnodes(parsize_*nvar);

  loccnt_ = 0;
  for( size_t n=0; n<parsize_; ++n )
  {
    if ( ! mask[n] )
    {
      sendnodes[loccnt_++] = glb_idx[n];
      sendnodes[loccnt_++] = part[n];
      sendnodes[loccnt_++] = remote_idx[n]-base;
    }
  }

  parallel::mpi::comm().allGather(loccnt_, glbcounts_.begin(), glbcounts_.end());

  glbcnt_ = std::accumulate(glbcounts_.begin(),glbcounts_.end(),0);

  glbdispls_[0]=0;
  for (size_t jproc = 1; jproc < nproc; ++jproc) // start at 1
  {
    glbdispls_[jproc]=glbcounts_[jproc-1]+glbdispls_[jproc-1];
  }
  std::vector<gidx_t> recvnodes(glbcnt_);

  parallel::mpi::comm().allGatherv(sendnodes.begin(), sendnodes.begin() + loccnt_,
                                recvnodes.data(), glbcounts_.data(), glbdispls_.data());

  // Load recvnodes in sorting structure
  size_t nb_recv_nodes = glbcnt_/nvar;
  std::vector<Node> node_sort(nb_recv_nodes);
  for( size_t n=0; n<nb_recv_nodes; ++n )
  {
    node_sort[n].g = recvnodes[n*nvar+0];
    node_sort[n].p = recvnodes[n*nvar+1];
    node_sort[n].i = recvnodes[n*nvar+2];
  }

  recvnodes.clear();

  // Sort on "g" member, and remove duplicates
  std::sort(node_sort.begin(), node_sort.end());
  node_sort.erase( std::unique( node_sort.begin(), node_sort.end() ), node_sort.end() );

  glbcounts_.assign(nproc,0);
  glbdispls_.assign(nproc,0);
  for( size_t n=0; n<node_sort.size(); ++n )
  {
    ++glbcounts_[node_sort[n].p] ;
  }
  glbdispls_[0]=0;

  for (size_t jproc = 1; jproc < nproc; ++jproc) // start at 1
  {
    glbdispls_[jproc]=glbcounts_[jproc-1]+glbdispls_[jproc-1];
  }
  glbcnt_ = std::accumulate(glbcounts_.begin(),glbcounts_.end(),0);
  loccnt_ = glbcounts_[myproc];


  glbmap_.clear(); glbmap_.resize(glbcnt_);
  locmap_.clear(); locmap_.resize(loccnt_);
  std::vector<int> idx(nproc,0);
  for( size_t n=0; n<node_sort.size(); ++n )
  {
    size_t jproc = node_sort[n].p;
    glbmap_[ glbdispls_[jproc]+idx[jproc] ] = n;

    if( jproc == myproc )
      locmap_[idx[jproc]] = node_sort[n].i;

    ++idx[jproc];
  }

  is_setup_ = true;
}



void GatherScatter::setup( const int part[],
                           const int remote_idx[], const int base,
                           const gidx_t glb_idx[],
                           const size_t parsize )
{
  std::vector<int> mask(parsize);
  IsGhostPoint is_ghost(part,remote_idx,base,parsize);
  for (size_t jj = 0; jj < parsize; ++jj)
  {
    mask[jj] =  is_ghost(jj) ? 1 : 0;
  }
  setup( part, remote_idx, base, glb_idx, mask.data(), parsize );
}

/////////////////////

GatherScatter* atlas__GatherScatter__new () {
  return new GatherScatter();
}

void atlas__GatherScatter__delete ( GatherScatter* This) {
  delete This;
}

void atlas__GatherScatter__setup32 ( GatherScatter* This,  int part[],
                                     int remote_idx[], int base,
                                     int glb_idx[],
                                     int parsize )
{
#if ATLAS_BITS_GLOBAL==32
  This->setup(part,remote_idx,base,glb_idx,parsize);
#else
  std::vector<gidx_t> glb_idx_convert(parsize);
  for(int j = 0; j < parsize; ++j)
  {
    glb_idx_convert[j] = glb_idx[j];
  }
  This->setup(part,remote_idx,base,glb_idx_convert.data(),parsize);
#endif
}

void atlas__GatherScatter__setup64 ( GatherScatter* This,  int part[],
                                     int remote_idx[], int base,
                                     long glb_idx[],
                                     int parsize )
{
#if ATLAS_BITS_GLOBAL==64
  This->setup(part,remote_idx,base,glb_idx,parsize);
#else
  std::vector<gidx_t> glb_idx_convert(parsize);
  for( size_t j=0; j<parsize; ++j )
  {
    glb_idx_convert[j] = glb_idx[j];
  }
  This->setup(part,remote_idx,base,glb_idx_convert.data(),parsize);
#endif
}

int atlas__GatherScatter__glb_dof ( GatherScatter* This )
{
  return This->glb_dof();
}

void atlas__GatherScatter__gather_int ( GatherScatter* This,
                                        int lfield[], int lvar_strides[], int lvar_extents[], int lvar_rank,
                                        int gfield[], int gvar_strides[], int gvar_extents[], int gvar_rank)
{
  std::vector<size_t> lvstrides(lvar_rank);
  std::vector<size_t> lvextents(lvar_rank);
  std::vector<size_t> gvstrides(gvar_rank);
  std::vector<size_t> gvextents(gvar_rank);
  for(int n = 0; n < lvar_rank; ++n) {
    lvstrides[n] = lvar_strides[n];
    lvextents[n] = lvar_extents[n];
  }
  for(int n = 0; n < gvar_rank; ++n) {
    gvstrides[n] = gvar_strides[n];
    gvextents[n] = gvar_extents[n];
  }
  This->gather(lfield,lvstrides.data(),lvextents.data(),lvar_rank,
               gfield,gvstrides.data(),gvextents.data(),gvar_rank);
}

void atlas__GatherScatter__gather_long ( GatherScatter* This,
                                         long lfield[], int lvar_strides[], int lvar_extents[], int lvar_rank,
                                         long gfield[], int gvar_strides[], int gvar_extents[], int gvar_rank)
{
  std::vector<size_t> lvstrides(lvar_rank);
  std::vector<size_t> lvextents(lvar_rank);
  std::vector<size_t> gvstrides(gvar_rank);
  std::vector<size_t> gvextents(gvar_rank);
  for(int n = 0; n < lvar_rank; ++n) {
    lvstrides[n] = lvar_strides[n];
    lvextents[n] = lvar_extents[n];
  }
  for(int n = 0; n < gvar_rank; ++n) {
    gvstrides[n] = gvar_strides[n];
    gvextents[n] = gvar_extents[n];
  }
  This->gather(lfield,lvstrides.data(),lvextents.data(),lvar_rank,
               gfield,gvstrides.data(),gvextents.data(),gvar_rank);}

void atlas__GatherScatter__gather_float ( GatherScatter* This,
                                          float lfield[], int lvar_strides[], int lvar_extents[], int lvar_rank,
                                          float gfield[], int gvar_strides[], int gvar_extents[], int gvar_rank)
{
  std::vector<size_t> lvstrides(lvar_rank);
  std::vector<size_t> lvextents(lvar_rank);
  std::vector<size_t> gvstrides(gvar_rank);
  std::vector<size_t> gvextents(gvar_rank);
  for(int n = 0; n < lvar_rank; ++n) {
    lvstrides[n] = lvar_strides[n];
    lvextents[n] = lvar_extents[n];
  }
  for(int n = 0; n < gvar_rank; ++n) {
    gvstrides[n] = gvar_strides[n];
    gvextents[n] = gvar_extents[n];
  }
  This->gather(lfield,lvstrides.data(),lvextents.data(),lvar_rank,
               gfield,gvstrides.data(),gvextents.data(),gvar_rank);}

void atlas__GatherScatter__gather_double ( GatherScatter* This,
                                           double lfield[], int lvar_strides[], int lvar_extents[], int lvar_rank,
                                           double gfield[], int gvar_strides[], int gvar_extents[], int gvar_rank)
{
  std::vector<size_t> lvstrides(lvar_rank);
  std::vector<size_t> lvextents(lvar_rank);
  std::vector<size_t> gvstrides(gvar_rank);
  std::vector<size_t> gvextents(gvar_rank);
  for(int n = 0; n < lvar_rank; ++n) {
    lvstrides[n] = lvar_strides[n];
    lvextents[n] = lvar_extents[n];
  }
  for(int n = 0; n < gvar_rank; ++n) {
    gvstrides[n] = gvar_strides[n];
    gvextents[n] = gvar_extents[n];
  }
  This->gather(lfield,lvstrides.data(),lvextents.data(),lvar_rank,
               gfield,gvstrides.data(),gvextents.data(),gvar_rank);}

void atlas__GatherScatter__scatter_int ( GatherScatter* This,
                                         int gfield[], int gvar_strides[], int gvar_extents[], int gvar_rank,
                                         int lfield[], int lvar_strides[], int lvar_extents[], int lvar_rank )
{
  std::vector<size_t> lvstrides(lvar_rank);
  std::vector<size_t> lvextents(lvar_rank);
  std::vector<size_t> gvstrides(gvar_rank);
  std::vector<size_t> gvextents(gvar_rank);
  for(int n = 0; n < lvar_rank; ++n) {
    lvstrides[n] = lvar_strides[n];
    lvextents[n] = lvar_extents[n];
  }
  for(int n = 0; n < gvar_rank; ++n) {
    gvstrides[n] = gvar_strides[n];
    gvextents[n] = gvar_extents[n];
  }
  This->scatter( gfield, gvstrides.data(), gvextents.data(), gvar_rank,
                 lfield, lvstrides.data(), lvextents.data(), lvar_rank );
}

void atlas__GatherScatter__scatter_long ( GatherScatter* This,
                                          long gfield[], int gvar_strides[], int gvar_extents[], int gvar_rank,
                                          long lfield[], int lvar_strides[], int lvar_extents[], int lvar_rank )
{
  std::vector<size_t> lvstrides(lvar_rank);
  std::vector<size_t> lvextents(lvar_rank);
  std::vector<size_t> gvstrides(gvar_rank);
  std::vector<size_t> gvextents(gvar_rank);
  for(int n = 0; n < lvar_rank; ++n) {
    lvstrides[n] = lvar_strides[n];
    lvextents[n] = lvar_extents[n];
  }
  for(int n = 0; n < gvar_rank; ++n) {
    gvstrides[n] = gvar_strides[n];
    gvextents[n] = gvar_extents[n];
  }
  This->scatter( gfield, gvstrides.data(), gvextents.data(), gvar_rank,
                 lfield, lvstrides.data(), lvextents.data(), lvar_rank );}

void atlas__GatherScatter__scatter_float ( GatherScatter* This,
                                           float gfield[], int gvar_strides[], int gvar_extents[], int gvar_rank,
                                           float lfield[], int lvar_strides[], int lvar_extents[], int lvar_rank )
{
  std::vector<size_t> lvstrides(lvar_rank);
  std::vector<size_t> lvextents(lvar_rank);
  std::vector<size_t> gvstrides(gvar_rank);
  std::vector<size_t> gvextents(gvar_rank);
  for(int n = 0; n < lvar_rank; ++n) {
    lvstrides[n] = lvar_strides[n];
    lvextents[n] = lvar_extents[n];
  }
  for(int n = 0; n < gvar_rank; ++n) {
    gvstrides[n] = gvar_strides[n];
    gvextents[n] = gvar_extents[n];
  }
  This->scatter( gfield, gvstrides.data(), gvextents.data(), gvar_rank,
                 lfield, lvstrides.data(), lvextents.data(), lvar_rank );}

void atlas__GatherScatter__scatter_double ( GatherScatter* This,
                                            double gfield[], int gvar_strides[], int gvar_extents[], int gvar_rank,
                                            double lfield[], int lvar_strides[], int lvar_extents[], int lvar_rank )
{
  std::vector<size_t> lvstrides(lvar_rank);
  std::vector<size_t> lvextents(lvar_rank);
  std::vector<size_t> gvstrides(gvar_rank);
  std::vector<size_t> gvextents(gvar_rank);
  for(int n = 0; n < lvar_rank; ++n) {
    lvstrides[n] = lvar_strides[n];
    lvextents[n] = lvar_extents[n];
  }
  for(int n = 0; n < gvar_rank; ++n) {
    gvstrides[n] = gvar_strides[n];
    gvextents[n] = gvar_extents[n];
  }
  This->scatter( gfield, gvstrides.data(), gvextents.data(), gvar_rank,
                 lfield, lvstrides.data(), lvextents.data(), lvar_rank );
}


/////////////////////

} // namespace parallel
} // namespace atlas
