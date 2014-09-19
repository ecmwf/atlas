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

#include "atlas/util/Array.h"
#include "atlas/util/ArrayView.h"
#include "atlas/util/Debug.h"
#include "atlas/mpl/GatherScatter.h"

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
    mypart_ = MPL::rank();
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
  Node() {}
  Node(int gid, int part, int idx)
  {
    g = gid;
    p = part;
    i = idx;
  }
  int g,p,i;
  bool operator < (const Node& other) const
  {
    return ( g < other.g );
  }
  bool operator == (const Node& other) const
  {
    return ( g == other.g );
  }
};


bool operator < (const int g, const Node& n)
{
  return ( g < n.g );
}

}

GatherScatter::GatherScatter() :
  is_setup_(false)
{
  myproc = MPL::rank();
  nproc  = MPL::size();
  root_   = 0;
}

void GatherScatter::setup( const int part[],
                           const int remote_idx[], const int base,
                           const int glb_idx[], const int max_glb_idx,
                           const int parsize )
{
  const bool include_ghost = false; // Warning: setting this to true might break scatter functionality
                                    //   This feature allows to check periodic halo values in output

  parsize_ = parsize;

  loccounts_.resize(nproc); loccounts_.assign(nproc,0);
  glbcounts_.resize(nproc); glbcounts_.assign(nproc,0);
  locdispls_.resize(nproc); locdispls_.assign(nproc,0);
  glbdispls_.resize(nproc); glbdispls_.assign(nproc,0);


  int maxgid = max_glb_idx;
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
    MPL_CHECK_RESULT( MPI_Allreduce(MPI_IN_PLACE,&maxgid,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD) );
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
  MPL_CHECK_RESULT( MPI_Gather( &loccnt_, 1, MPI_INT,
                     glbcounts_.data(), 1, MPI_INT,
                     root_, MPI_COMM_WORLD ) );
  glbcnt_ = std::accumulate(glbcounts_.begin(),glbcounts_.end(),0);

  glbdispls_[0]=0;
  for (int jproc=1; jproc<nproc; ++jproc) // start at 1
  {
    glbdispls_[jproc]=glbcounts_[jproc-1]+glbdispls_[jproc-1];
  }
  std::vector<int> recvnodes(glbcnt_);
  MPL_CHECK_RESULT( MPI_Gatherv( sendnodes.data(), loccnt_, MPI_INT,
                                 recvnodes.data(), glbcounts_.data(), glbdispls_.data(), MPI_INT,
                      root_, MPI_COMM_WORLD) );

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
  MPL_CHECK_RESULT( MPI_Scatter( glbcounts_.data(), 1, MPI_INT,
                                 &loccnt_,     1, MPI_INT,
                                 root_, MPI_COMM_WORLD) );

//  DEBUG_VAR_SYNC(loccnt_);
  locmap_.resize(loccnt_);

  MPL_CHECK_RESULT( MPI_Scatterv( needed.data(), glbcounts_.data(), glbdispls_.data(),
                                  MPI_INT, locmap_.data(), loccnt_,
                                  MPI_INT, root_, MPI_COMM_WORLD ) );


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

void atlas__GatherScatter__setup ( GatherScatter* This,  int part[],
                                   int remote_idx[], int base,
                                   int glb_idx[], int max_glb_idx,
                                   int parsize )
{
  This->setup(part,remote_idx,base,glb_idx,max_glb_idx,parsize);
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
