#include <stdexcept>
#include <numeric>      // std::accumulate
#include <iostream>
#include <sstream>
#include "Comm.hpp"

namespace ecmwf {

HaloExchange::HaloExchange()
{
  int ierr;
  ierr = MPI_Comm_rank( MPI_COMM_WORLD, &myproc );
  ierr = MPI_Comm_size( MPI_COMM_WORLD, &nproc );
}

void HaloExchange::setup( const int proc[], 
                          const int glb_idx[], 
                          const std::vector<int>& bounds, 
                          int par_bound )
{

  int ierr;


//  bounds_.resize(bounds.size());
//  for(int i=0; i<bounds.size(); ++i)
//    bounds_[i] = bounds[i];

  bounds_ = bounds;
  par_bound_ = par_bound;


  sync_sendcounts_.resize(nproc,0);
  sync_recvcounts_.resize(nproc,0);
  sync_senddispls_.resize(nproc,0);
  sync_recvdispls_.resize(nproc,0);

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

  ierr = MPI_Allreduce( MPI_IN_PLACE, &max_glb_idx, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD );

  std::vector<int> map_glb_to_loc(max_glb_idx+1,-1);
  for (int jj=0; jj<nb_nodes; ++jj)
    map_glb_to_loc[glb_idx[jj]] = jj;

  // std::cout << myproc << ":  map_glb_to_loc = ";
  // for (int i=0; i<max_glb_idx+1; ++i)
  //   std::cout << map_glb_to_loc[i] << " ";
  // std::cout << std::endl;

  /*
    Find the amount of nodes this proc has to receive from each other proc
  */

  for (int jj=0; jj<nb_nodes; ++jj)
  {
    if (proc[jj] != myproc)
      ++sync_recvcounts_[proc[jj]];
  }
  sync_recvcnt_ = std::accumulate(sync_recvcounts_.begin(),sync_recvcounts_.end(),0);
  // std::cout << myproc << ":  recvcnt = " << sync_recvcnt_ << std::endl;


  /*
    Find the amount of nodes this proc has to send to each other proc
  */

  ierr = MPI_Alltoall( &sync_recvcounts_[0], 1, MPI_INT, &sync_sendcounts_[0], 1, MPI_INT, MPI_COMM_WORLD );
  sync_sendcnt_ = std::accumulate(sync_sendcounts_.begin(),sync_sendcounts_.end(),0);

  sync_recvdispls_[0]=0;
  sync_senddispls_[0]=0;
  for (int jproc=1; jproc<nproc; ++jproc) // start at 1
  {
    sync_recvdispls_[jproc]=sync_recvcounts_[jproc-1]+sync_recvdispls_[jproc-1];
    sync_senddispls_[jproc]=sync_sendcounts_[jproc-1]+sync_senddispls_[jproc-1];
  }
  /*
    Fill vector "send_requests" with global indices of nodes needed, but are on other procs
    We can also fill in the vector "sync_recvmap_" which holds local indices of requested nodes
  */

  std::vector<int> send_requests(sync_recvcnt_);

  sync_recvmap_.resize(sync_recvcnt_);
  std::vector<int> cnt(nproc,0);
  for (int jj=0; jj<nb_nodes; ++jj)
  {
    if (proc[jj] != myproc)
    {
      const int req_idx = sync_recvdispls_[proc[jj]] + cnt[proc[jj]];
      send_requests[req_idx] = glb_idx[jj];
      sync_recvmap_[req_idx] = jj;
      ++cnt[proc[jj]];
    }
  }

  /*
    Fill vector "recv_requests" with global indices that are needed by other procs
  */

  std::vector<int> recv_requests(sync_sendcnt_);

  ierr = MPI_Alltoallv( &send_requests[0], &sync_recvcounts_[0], &sync_recvdispls_[0], MPI_INT,
                        &recv_requests[0], &sync_sendcounts_[0], &sync_senddispls_[0], MPI_INT,
                        MPI_COMM_WORLD );

  /*
    What needs to be sent to other procs can be found by a map from global to local indices
  */
  sync_sendmap_.resize(sync_sendcnt_);
  for( int jj=0; jj<sync_sendcnt_; ++jj )
    sync_sendmap_[jj] = map_glb_to_loc[ recv_requests[jj] ];

  // Packet size
  packet_size_ = 1;
  const int nb_bounds = bounds_.size();
  for (int b=0; b<nb_bounds; ++b)
  {
    if ( b != par_bound_ && bounds_[b] >= 0)
      packet_size_ *= bounds_[b];
  }

  // std::cout << myproc << "  :  sync_sendmap_  = ";
  // for( int i=0; i< sync_sendmap_.size(); ++i)
  //   std::cout << sync_sendmap_[i] << " ";
  // std::cout << std::endl;

  // std::cout << myproc << "  :  sync_recvmap_  = ";
  // for( int i=0; i< sync_recvmap_.size(); ++i)
  //   std::cout << sync_recvmap_[i] << " ";
  // std::cout << std::endl;
}



template<>
inline void HaloExchange::create_mappings_impl<2,0>( 
    std::vector<int>& send_map, 
    std::vector<int>& recv_map,
    int nb_vars) const
{
  const int nb_nodes = bounds_[0];
  int send_idx(0);
  for (int jnode=0; jnode<sync_sendcnt_; ++jnode)
  {
    const int inode = sync_sendmap_[jnode];
    for (int jvar=0; jvar<nb_vars; ++jvar)
    {
      send_map[send_idx++] = index( inode,jvar,   nb_nodes,nb_vars);
    }
  }
  int recv_idx(0);
  for (int jnode=0; jnode<sync_recvcnt_; ++jnode)
  {
    const int inode = sync_recvmap_[jnode];
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
  for (int jnode=0; jnode<sync_sendcnt_; ++jnode)
  {
    const int inode = sync_sendmap_[jnode];
    for (int jvar=0; jvar<nb_vars; ++jvar)
    {
      send_map[send_idx++] = index( jvar, inode,  nb_vars, nb_nodes);
    }
  }
  int recv_idx(0);
  for (int jnode=0; jnode<sync_recvcnt_; ++jnode)
  {
    const int inode = sync_recvmap_[jnode];
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
  for (int n=0; n<sync_sendcnt_; ++n)
  {
    const int jnode = sync_sendmap_[n];
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

  for (int n=0; n<sync_recvcnt_; ++n)
  {
    const int jnode = sync_recvmap_[n];
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


HaloExchange* ecmwf__HaloExchange__new () { 
  return new HaloExchange(); 
}

void ecmwf__HaloExchange__delete (HaloExchange* This) {
  delete This;
}

void ecmwf__HaloExchange__setup (HaloExchange* This, int proc[], int glb_idx[], int bounds[], int nb_bounds, int par_bound)
{
  std::vector<int> bounds_vec(bounds,bounds+nb_bounds);
  This->setup(proc,glb_idx,bounds_vec,par_bound);
}

void ecmwf__HaloExchange__execute_int (HaloExchange* This, int field[], int nb_vars ) { 
  This->execute(field,nb_vars);
}

void ecmwf__HaloExchange__execute_float (HaloExchange* This, float field[], int nb_vars ) { 
  This->execute(field,nb_vars);
}

void ecmwf__HaloExchange__execute_double (HaloExchange* This, double field[], int nb_vars ) { 
  This->execute(field,nb_vars);
}

/////////////////////

}
