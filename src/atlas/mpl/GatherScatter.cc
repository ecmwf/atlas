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

#include "eckit/log/Log.h"
#include "atlas/util/Array.h"
#include "atlas/util/ArrayView.h"
#include "atlas/util/Debug.h"
#include "atlas/mpl/GatherScatter.h"
#include "atlas/util/Checksum.h"

using eckit::Log;

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
  is_setup_(false)
{
  myproc = eckit::mpi::rank();
  nproc  = eckit::mpi::size();
  root_   = 0;
}



void GatherScatter::setup( const int part[],
                           const int remote_idx[], const int base,
                           const gidx_t glb_idx[], const int mask[], const int parsize )
{
  parsize_ = parsize;

  loccounts_.resize(nproc); loccounts_.assign(nproc,0);
  glbcounts_.resize(nproc); glbcounts_.assign(nproc,0);
  locdispls_.resize(nproc); locdispls_.assign(nproc,0);
  glbdispls_.resize(nproc); glbdispls_.assign(nproc,0);

  const int nvar = 3;

  std::vector<gidx_t> sendnodes(parsize_*nvar);

  loccnt_ = 0;
  for( int n=0; n<parsize_; ++n )
  {
    if ( ! mask[n] )
    {
      sendnodes[loccnt_++] = glb_idx[n];
      sendnodes[loccnt_++] = part[n];
      sendnodes[loccnt_++] = remote_idx[n]-base;
    }
  }

  ECKIT_MPI_CHECK_RESULT( MPI_Gather( &loccnt_, 1, MPI_INT,
                                glbcounts_.data(), 1, MPI_INT,
                                root_, eckit::mpi::comm() ) );
  glbcnt_ = std::accumulate(glbcounts_.begin(),glbcounts_.end(),0);

  glbdispls_[0]=0;
  for (int jproc=1; jproc<nproc; ++jproc) // start at 1
  {
    glbdispls_[jproc]=glbcounts_[jproc-1]+glbdispls_[jproc-1];
  }
  std::vector<gidx_t> recvnodes(glbcnt_);
  ECKIT_MPI_CHECK_RESULT( MPI_Gatherv( sendnodes.data(), loccnt_, eckit::mpi::datatype<gidx_t>(),
                                 recvnodes.data(), glbcounts_.data(), glbdispls_.data(), eckit::mpi::datatype<gidx_t>(),
                                 root_, eckit::mpi::comm()) );

  // Load recvnodes in sorting structure
  int nb_recv_nodes = glbcnt_/nvar;
  std::vector<Node> node_sort(nb_recv_nodes);
  for( int n=0; n<nb_recv_nodes; ++n )
  {
    node_sort[n].g = recvnodes[n*3+0];
    node_sort[n].p = recvnodes[n*3+1];
    node_sort[n].i = recvnodes[n*3+2];
  }

  recvnodes.clear();

  // Sort on "g" member, and remove duplicates
  std::sort(node_sort.begin(), node_sort.end());
  node_sort.erase( std::unique( node_sort.begin(), node_sort.end() ), node_sort.end() );

  glbcounts_.assign(nproc,0);
  glbdispls_.assign(nproc,0);
  for( int n=0; n<node_sort.size(); ++n )
  {
    ++glbcounts_[node_sort[n].p] ;
  }
  glbdispls_[0]=0;

  for (int jproc=1; jproc<nproc; ++jproc) // start at 1
  {
    glbdispls_[jproc]=glbcounts_[jproc-1]+glbdispls_[jproc-1];
  }
  glbcnt_ = std::accumulate(glbcounts_.begin(),glbcounts_.end(),0);


  glbmap_.clear(); glbmap_.resize(glbcnt_);
  std::vector<int> needed(glbcnt_);
  std::vector<int> idx(nproc,0);
  for( int n=0; n<node_sort.size(); ++n )
  {
    int jproc = node_sort[n].p;
    needed [ glbdispls_[jproc]+idx[jproc] ] = node_sort[n].i; // index on sending proc
    glbmap_[ glbdispls_[jproc]+idx[jproc] ] = n;
    ++idx[jproc];
  }

  // Get loccnt_
  ECKIT_MPI_CHECK_RESULT( MPI_Scatter( glbcounts_.data(), 1, MPI_INT,
                                 &loccnt_,          1, MPI_INT,
                                 root_, eckit::mpi::comm()) );

  locmap_.resize(loccnt_);

  ECKIT_MPI_CHECK_RESULT( MPI_Scatterv( needed.data(), glbcounts_.data(), glbdispls_.data(),
                                  MPI_INT, locmap_.data(), loccnt_,
                                  MPI_INT, root_, eckit::mpi::comm() ) );
  is_setup_ = true;
}



void GatherScatter::setup( const int part[],
                           const int remote_idx[], const int base,
                           const gidx_t glb_idx[], const gidx_t max_glb_idx,
                           const int parsize, const bool include_ghost )
{
  parsize_ = parsize;

  loccounts_.resize(nproc); loccounts_.assign(nproc,0);
  glbcounts_.resize(nproc); glbcounts_.assign(nproc,0);
  locdispls_.resize(nproc); locdispls_.assign(nproc,0);
  glbdispls_.resize(nproc); glbdispls_.assign(nproc,0);


  gidx_t maxgid = max_glb_idx;
  if(max_glb_idx<0)
  {
    IsGhostPoint is_ghost(part,remote_idx,base,parsize_);
    maxgid = -1;
    for (int jj=0; jj<parsize_; ++jj)
    {
      if ( !is_ghost(jj) || include_ghost )
      {
        maxgid = std::max(maxgid,glb_idx[jj]);
      }
    }
    ECKIT_MPI_CHECK_RESULT( MPI_Allreduce(MPI_IN_PLACE,&maxgid,1,eckit::mpi::datatype<gidx_t>(),MPI_MAX,eckit::mpi::comm()) );
  }

  Array<int> sendnodes(parsize_,3);
  ArrayView<int,2> nodes(sendnodes);

  if( include_ghost )
  {
    for( int n=0; n<parsize_; ++n )
    {
      nodes(n,0) = glb_idx[n];
      nodes(n,1) = myproc;
      nodes(n,2) = n;
    }
  }
  else
  {
    for( int n=0; n<parsize_; ++n )
    {
      nodes(n,0) = glb_idx[n];
      nodes(n,1) = part[n];
      nodes(n,2) = remote_idx[n]-base;
    }

  }

  loccnt_ = nodes.total_size();
  ECKIT_MPI_CHECK_RESULT( MPI_Gather( &loccnt_, 1, MPI_INT,
                     glbcounts_.data(), 1, MPI_INT,
                     root_, eckit::mpi::comm() ) );
  glbcnt_ = std::accumulate(glbcounts_.begin(),glbcounts_.end(),0);

  glbdispls_[0]=0;
  for (int jproc=1; jproc<nproc; ++jproc) // start at 1
  {
    glbdispls_[jproc]=glbcounts_[jproc-1]+glbdispls_[jproc-1];
  }
  std::vector<int> recvnodes(glbcnt_);
  ECKIT_MPI_CHECK_RESULT( MPI_Gatherv( sendnodes.data(), loccnt_, MPI_INT,
                                 recvnodes.data(), glbcounts_.data(), glbdispls_.data(), MPI_INT,
                      root_, eckit::mpi::comm()) );

  // Load recvnodes in sorting structure
  int nb_recv_nodes = glbcnt_/3;
  std::vector<Node> node_sort(nb_recv_nodes);
  nodes = ArrayView<int,2> (recvnodes.data(),make_shape(nb_recv_nodes,3));
  for( int n=0; n<nb_recv_nodes; ++n )
  {
    node_sort[n].g = nodes(n,0);
    node_sort[n].p = nodes(n,1);
    node_sort[n].i = nodes(n,2);
  }

  recvnodes.clear();

  // Sort on "g" member, remove nodes with g larger than max_glb_idx, and remove duplicates
  std::sort(node_sort.begin(), node_sort.end());
  node_sort.erase( std::upper_bound ( node_sort.begin(), node_sort.end(), maxgid ), node_sort.end() );
  node_sort.erase( std::unique( node_sort.begin(), node_sort.end() ), node_sort.end() );


//    if ( myproc == root_ )
//    {
//      std::vector<gidx_t> gnodes(node_sort.size());
//      for( int i=0; i<node_sort.size(); ++i)
//      {
//        gnodes[i] = node_sort[i].g;
////        if(!include_ghost)
////        {
////        //std::cout << "gnodes["<<i<<"] = " << gnodes[i] << std::endl;
////        ASSERT( gnodes[i] == i+1 );
////        }
//      }
//      //eckit::Log::info() << "checksum glb_idx[0:"<<node_sort.size()<<"] = " << checksum(gnodes.data(),gnodes.size()) << std::endl;
//    }

//  if ( myproc == root )
//  {
//    std::cout << myproc << "  :  node_sort  = ";
//    for( int i=0; i< node_sort.size(); ++i)
//      std::cout << node_sort[i].g << " ";
//    std::cout << std::endl;
//  }

  // Assemble list to ask needed
  glbcounts_.assign(nproc,0);
  glbdispls_.assign(nproc,0);
  for( int n=0; n<node_sort.size(); ++n )
  {
    ++glbcounts_[node_sort[n].p] ;
  }
  glbdispls_[0]=0;

  for (int jproc=1; jproc<nproc; ++jproc) // start at 1
  {
    glbdispls_[jproc]=glbcounts_[jproc-1]+glbdispls_[jproc-1];
  }
  glbcnt_ = std::accumulate(glbcounts_.begin(),glbcounts_.end(),0);


  glbmap_.clear(); glbmap_.resize(glbcnt_);
  std::vector<int> needed(glbcnt_);
  std::vector<int> idx(nproc,0);
  for( int n=0; n<node_sort.size(); ++n )
  {
    int jproc = node_sort[n].p;
    needed[ glbdispls_[jproc]+idx[jproc] ] = node_sort[n].i; // index on sending proc
    glbmap_[glbdispls_[jproc]+idx[jproc]] = n;
    ++idx[jproc];
  }

//  if ( myproc == root ) std::cout << "needed = ";
//  for( int n=0; n<glbcnt_; ++n )
//  {
//    if( myproc == root ) std::cout << needed[n] << " " ;
//  }
//  if( myproc == root ) std::cout << std::endl;




//  std::cout << myproc << "  :  glbcounts_  = ";
//  for( int i=0; i< nproc; ++i)
//    std::cout << glbcounts_[i] << " ";
//  std::cout << std::endl;
//  std::cout << myproc << "  :  glbdispls_  = ";
//  for( int i=0; i< nproc; ++i)
//    std::cout << glbdispls_[i] << " ";
//  std::cout << std::endl;

  // Get loccnt_
  ECKIT_MPI_CHECK_RESULT( MPI_Scatter( glbcounts_.data(), 1, MPI_INT,
                                 &loccnt_,     1, MPI_INT,
                                 root_, eckit::mpi::comm()) );

//  DEBUG_VAR_SYNC(loccnt_);
  locmap_.resize(loccnt_);

  ECKIT_MPI_CHECK_RESULT( MPI_Scatterv( needed.data(), glbcounts_.data(), glbdispls_.data(),
                                  MPI_INT, locmap_.data(), loccnt_,
                                  MPI_INT, root_, eckit::mpi::comm() ) );


//   std::cout << myproc << "  :  locmap_  = ";
//   for( int i=0; i< loccnt_; ++i)
//     std::cout << locmap_[i] << " ";
//   std::cout << std::endl;


//  for( int i=0; i<glbcnt_; ++i )
//  {
//    ASSERT( node_sort[i].g == i+1 );
//  }


  is_setup_ = true;
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
                                     int glb_idx[], int max_glb_idx,
                                     int parsize )
{
#if ATLAS_BITS_GLOBAL==32
  This->setup(part,remote_idx,base,glb_idx,max_glb_idx,parsize);
#else
  std::vector<gidx_t> glb_idx_convert(parsize);
  for( int j=0; j<parsize; ++j )
  {
    glb_idx_convert[j] = glb_idx[j];
  }
  This->setup(part,remote_idx,base,glb_idx_convert.data(),max_glb_idx,parsize);
#endif
}

void atlas__GatherScatter__setup64 ( GatherScatter* This,  int part[],
                                     int remote_idx[], int base,
                                     long glb_idx[], long max_glb_idx,
                                     int parsize )
{
#if ATLAS_BITS_GLOBAL==64
  This->setup(part,remote_idx,base,glb_idx,max_glb_idx,parsize);
#else
  std::vector<gidx_t> glb_idx_convert(parsize);
  for( int j=0; j<parsize; ++j )
  {
    glb_idx_convert[j] = glb_idx[j];
  }
  This->setup(part,remote_idx,base,glb_idx_convert.data(),max_glb_idx,parsize);
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
  This->gather(lfield,lvar_strides,lvar_extents,lvar_rank,
               gfield,gvar_strides,gvar_extents,gvar_rank);
}

void atlas__GatherScatter__gather_long ( GatherScatter* This,
                                         long lfield[], int lvar_strides[], int lvar_extents[], int lvar_rank,
                                         long gfield[], int gvar_strides[], int gvar_extents[], int gvar_rank)
{
  This->gather(lfield,lvar_strides,lvar_extents,lvar_rank,
               gfield,gvar_strides,gvar_extents,gvar_rank);
}

void atlas__GatherScatter__gather_float ( GatherScatter* This,
                                          float lfield[], int lvar_strides[], int lvar_extents[], int lvar_rank,
                                          float gfield[], int gvar_strides[], int gvar_extents[], int gvar_rank)
{
  This->gather(lfield,lvar_strides,lvar_extents,lvar_rank,
               gfield,gvar_strides,gvar_extents,gvar_rank);
}

void atlas__GatherScatter__gather_double ( GatherScatter* This,
                                           double lfield[], int lvar_strides[], int lvar_extents[], int lvar_rank,
                                           double gfield[], int gvar_strides[], int gvar_extents[], int gvar_rank)
{
  This->gather(lfield,lvar_strides,lvar_extents,lvar_rank,
               gfield,gvar_strides,gvar_extents,gvar_rank);
}

void atlas__GatherScatter__scatter_int ( GatherScatter* This,
                                         int gfield[], int gvar_strides[], int gvar_extents[], int gvar_rank,
                                         int lfield[], int lvar_strides[], int lvar_extents[], int lvar_rank )
{
  This->scatter( gfield, gvar_strides, gvar_extents, gvar_rank,
                 lfield, lvar_strides, lvar_extents, lvar_rank );
}

void atlas__GatherScatter__scatter_long ( GatherScatter* This,
                                          long gfield[], int gvar_strides[], int gvar_extents[], int gvar_rank,
                                          long lfield[], int lvar_strides[], int lvar_extents[], int lvar_rank )
{
  This->scatter( gfield, gvar_strides, gvar_extents, gvar_rank,
                 lfield, lvar_strides, lvar_extents, lvar_rank );
}

void atlas__GatherScatter__scatter_float ( GatherScatter* This,
                                           float gfield[], int gvar_strides[], int gvar_extents[], int gvar_rank,
                                           float lfield[], int lvar_strides[], int lvar_extents[], int lvar_rank )
{
  This->scatter( gfield, gvar_strides, gvar_extents, gvar_rank,
                 lfield, lvar_strides, lvar_extents, lvar_rank );
}

void atlas__GatherScatter__scatter_double ( GatherScatter* This,
                                            double gfield[], int gvar_strides[], int gvar_extents[], int gvar_rank,
                                            double lfield[], int lvar_strides[], int lvar_extents[], int lvar_rank )
{
  This->scatter( gfield, gvar_strides, gvar_extents, gvar_rank,
                 lfield, lvar_strides, lvar_extents, lvar_rank );
}


/////////////////////
} // namespace mpl
} // namespace atlas
