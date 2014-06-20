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

#include "atlas/util/Array.hpp"
#include "atlas/util/ArrayView.hpp"
#include "atlas/util/Debug.hpp"
#include "atlas/mpl/Gather.hpp"

namespace atlas {

namespace {

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
}

Gather::Gather() :
  is_setup_(false)
{
  myproc = MPL::rank();
  nproc  = MPL::size();
  root   = 0;
}

void Gather::setup( const int part[],
                    const int remote_idx[], const int base,
                    const int glb_idx[], const int max_glb_idx,
                    const int parsize )
{
  parsize_ = parsize;

  sendcounts_.resize(nproc,0);
  recvcounts_.resize(nproc,0);
  senddispls_.resize(nproc,0);
  recvdispls_.resize(nproc,0);

  Array<int> sendnodes(parsize_,3);
  ArrayView<int,2> nodes(sendnodes);
  for( int n=0; n<parsize_; ++n )
  {
    nodes(n,0) = glb_idx[n];
    nodes(n,1) = part[n];
    nodes(n,2) = remote_idx[n]-base;
  }

  sendcnt_ = nodes.total_size();
  MPL_CHECK_RESULT( MPI_Gather( &sendcnt_, 1, MPI_INT,
                     recvcounts_.data(), 1, MPI_INT,
                     root, MPI_COMM_WORLD ) );
  recvcnt_ = std::accumulate(recvcounts_.begin(),recvcounts_.end(),0);

  recvdispls_[0]=0;
  for (int jproc=1; jproc<nproc; ++jproc) // start at 1
  {
    recvdispls_[jproc]=recvcounts_[jproc-1]+recvdispls_[jproc-1];
  }
  std::vector<int> recvnodes(recvcnt_);
  MPL_CHECK_RESULT( MPI_Gatherv( sendnodes.data(), sendcnt_, MPI_INT,
                                 recvnodes.data(), recvcounts_.data(), recvdispls_.data(), MPI_INT,
                      root, MPI_COMM_WORLD) );

  // Load recvnodes in sorting structure
  int nb_recv_nodes = recvcnt_/3;
  std::vector<Node> node_sort(nb_recv_nodes);
  nodes = ArrayView<int,2> (recvnodes.data(),Extents(nb_recv_nodes,3));
  for( int n=0; n<nb_recv_nodes; ++n )
  {
    node_sort[n].g = nodes(n,0);
    node_sort[n].p = nodes(n,1);
    node_sort[n].i = nodes(n,2);
  }

  recvnodes.clear();

  // Sort on "g" member, remove nodes with g larger than max_glb_idx, and remove duplicates
  std::sort(node_sort.begin(), node_sort.end());
  node_sort.erase( std::upper_bound ( node_sort.begin(), node_sort.end(), max_glb_idx ), node_sort.end() );
  node_sort.erase( std::unique( node_sort.begin(), node_sort.end() ), node_sort.end() );

//  if ( myproc == root )
//  {
//    std::cout << myproc << "  :  node_sort  = ";
//    for( int i=0; i< node_sort.size(); ++i)
//      std::cout << node_sort[i].g << " ";
//    std::cout << std::endl;
//  }

  // Assemble list to ask needed
  recvcounts_.assign(nproc,0);
  recvdispls_.assign(nproc,0);
  for( int n=0; n<node_sort.size(); ++n )
  {
    ++recvcounts_[node_sort[n].p] ;
  }
  recvdispls_[0]=0;
  for (int jproc=1; jproc<nproc; ++jproc) // start at 1
  {
    recvdispls_[jproc]=recvcounts_[jproc-1]+recvdispls_[jproc-1];
  }
  recvcnt_ = std::accumulate(recvcounts_.begin(),recvcounts_.end(),0);

  std::vector<int> needed(recvcnt_);
  std::vector<int> idx(nproc,0);
  for( int n=0; n<node_sort.size(); ++n )
  {
    int jproc = node_sort[n].p;
    needed[ recvdispls_[jproc]+idx[jproc] ] = node_sort[n].i; // index on sending proc
    ++idx[jproc];
  }

//  if ( myproc == root ) std::cout << "needed = ";
//  for( int n=0; n<recvcnt_; ++n )
//  {
//    if( myproc == root ) std::cout << needed[n] << " " ;
//  }
//  if( myproc == root ) std::cout << std::endl;




//  std::cout << myproc << "  :  recvcounts_  = ";
//  for( int i=0; i< nproc; ++i)
//    std::cout << recvcounts_[i] << " ";
//  std::cout << std::endl;
//  std::cout << myproc << "  :  recvdispls_  = ";
//  for( int i=0; i< nproc; ++i)
//    std::cout << recvdispls_[i] << " ";
//  std::cout << std::endl;

  // Get sendcnt_
  MPL_CHECK_RESULT( MPI_Scatter( recvcounts_.data(), 1, MPI_INT,
                                 &sendcnt_,     1, MPI_INT,
                                 root, MPI_COMM_WORLD) );

//  DEBUG_VAR_SYNC(sendcnt_);
  sendmap_.resize(sendcnt_);

  MPL_CHECK_RESULT( MPI_Scatterv( needed.data(), recvcounts_.data(), recvdispls_.data(),
                                  MPI_INT, sendmap_.data(), sendcnt_,
                                  MPI_INT, root, MPI_COMM_WORLD ) );


//   std::cout << myproc << "  :  sendmap_  = ";
//   for( int i=0; i< sendcnt_; ++i)
//     std::cout << sendmap_[i] << " ";
//   std::cout << std::endl;


  is_setup_ = true;
}


/////////////////////


Gather* atlas__Gather__new () {
  return new Gather();
}

void atlas__Gather__delete (Gather* This) {
  delete This;
}

void atlas__Gather__setup (Gather* This,  int part[],
                           int remote_idx[], int base,
                           int glb_idx[], int max_glb_idx,
                           int parsize )
{
  This->setup(part,remote_idx,base,glb_idx,max_glb_idx,parsize);
}

void atlas__Gather__execute_strided_int (Gather* This,
                                         int lfield[], int lvar_strides[], int lvar_extents[], int lvar_rank,
                                         int gfield[], int gvar_strides[], int gvar_extents[], int gvar_rank)
{
  This->execute(lfield,lvar_strides,lvar_extents,lvar_rank,
                gfield,gvar_strides,gvar_extents,gvar_rank);
}

void atlas__Gather__execute_strided_float (Gather* This,
                                           float lfield[], int lvar_strides[], int lvar_extents[], int lvar_rank,
                                           float gfield[], int gvar_strides[], int gvar_extents[], int gvar_rank)
{
  This->execute(lfield,lvar_strides,lvar_extents,lvar_rank,
                gfield,gvar_strides,gvar_extents,gvar_rank);
}

void atlas__Gather__execute_strided_double (Gather* This,
                                            double lfield[], int lvar_strides[], int lvar_extents[], int lvar_rank,
                                            double gfield[], int gvar_strides[], int gvar_extents[], int gvar_rank)
{
  This->execute(lfield,lvar_strides,lvar_extents,lvar_rank,
                gfield,gvar_strides,gvar_extents,gvar_rank);
}

void atlas__Gather__execute_int (Gather* This, int locfield[], int glbfield[], int nb_vars ) {
  This->execute(locfield,glbfield,nb_vars);
}

void atlas__Gather__execute_float (Gather* This, float locfield[], float glbfield[], int nb_vars ) {
  This->execute(locfield,glbfield,nb_vars);
}

void atlas__Gather__execute_double (Gather* This, double locfield[], double glbfield[], int nb_vars ) {
  This->execute(locfield,glbfield,nb_vars);
}

/////////////////////

}
