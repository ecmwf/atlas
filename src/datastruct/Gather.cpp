#include <stdexcept>
#include <numeric>      // std::accumulate
#include <iostream>
#include <sstream>
#include "Gather.hpp"

namespace ecmwf {

Gather::Gather() :
  is_setup_(false)
{
  myproc = MPL::rank();
  nproc  = MPL::size();
  root   = 0;
}

void Gather::setup(const int proc[],
                   const int glb_idx[],
                   const int master_glb_idx[],
                   const std::vector<int>& bounds,
                   int par_bound )
{

  int ierr;

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

//  int max_glb_idx = -1;
//  for (int jj=0; jj<nb_nodes; ++jj)
//  {
//    max_glb_idx = std::max( max_glb_idx, glb_idx[jj] );
//  }

//  ierr = MPI_Allreduce( MPI_IN_PLACE, &max_glb_idx, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD );

  /*
    Find the amount of nodes this proc has to send
  */


  sendcnt_ = 0;
  for (int jj=0; jj<nb_nodes; ++jj)
  {
    if (proc[jj] == myproc)
      ++sendcnt_;
  }

  ierr = MPI_Gather( &sendcnt_, 1, MPI_INT,
                     recvcounts_.data(), 1, MPI_INT,
                     root, MPI_COMM_WORLD );

  recvcnt_ = std::accumulate(recvcounts_.begin(),recvcounts_.end(),0);

  recvdispls_[0]=0;
  for (int jproc=1; jproc<nproc; ++jproc) // start at 1
  {
    recvdispls_[jproc]=recvcounts_[jproc-1]+recvdispls_[jproc-1];
  }

  /*
    What needs to be sent to other procs can be found by a map from global to local indices
  */
  std::vector<int> send_recv_position(sendcnt_);
  sendmap_.resize(sendcnt_);
  recvmap_.resize(recvcnt_);

  int idx_send(0);
  for (int jnode=0; jnode<nb_nodes; ++jnode)
  {
    if (proc[jnode] == myproc)
    {
      sendmap_[idx_send] = jnode;
      send_recv_position[idx_send] = glb_idx[jnode]-1;
      ++idx_send;
    }
  }

  ierr = MPI_Gatherv( send_recv_position.data(), sendcnt_, MPI_INT,
                      recvmap_.data(), recvcounts_.data(), recvdispls_.data(), MPI_INT,
                      root, MPI_COMM_WORLD);


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


/////////////////////


Gather* ecmwf__Gather__new () {
  return new Gather();
}

void ecmwf__Gather__delete (Gather* This) {
  delete This;
}

void ecmwf__Gather__setup (Gather* This, int proc[], int glb_idx[], int master_glb_idx[], int bounds[], int nb_bounds, int par_bound)
{
  std::vector<int> bounds_vec(bounds,bounds+nb_bounds);
  This->setup(proc,glb_idx,master_glb_idx,bounds_vec,par_bound);
}

void ecmwf__Gather__execute_int (Gather* This, int locfield[], int glbfield[], int nb_vars ) {
  This->execute(locfield,glbfield,nb_vars);
}

void ecmwf__Gather__execute_float (Gather* This, float locfield[], float glbfield[], int nb_vars ) {
  This->execute(locfield,glbfield,nb_vars);
}

void ecmwf__Gather__execute_double (Gather* This, double locfield[], double glbfield[], int nb_vars ) {
  This->execute(locfield,glbfield,nb_vars);
}

/////////////////////

}
