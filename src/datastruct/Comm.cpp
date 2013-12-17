#include <stdexcept>
#include <numeric>      // std::accumulate
#include <iostream>
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
  int max_glb_idx;
  int ierr;
  int nb_nodes;

  bounds_ = bounds;
  par_bound_ = par_bound;

  sync_sendcounts_.resize(nproc,0);
  sync_recvcounts_.resize(nproc,0);
  sync_senddispls_.resize(nproc,0);
  sync_recvdispls_.resize(nproc,0);

  nb_nodes = bounds_[par_bound_];

  for (int jj=0; jj<nb_nodes; ++jj)
  {
    max_glb_idx = std::max( max_glb_idx, glb_idx[jj] );
  }

  ierr = MPI_Allreduce( MPI_IN_PLACE, &max_glb_idx, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD );

  std::vector<int> map_glb_to_loc(max_glb_idx+1); 
  for (int jj=0; jj<nb_nodes; ++jj)
  {
    if (proc[jj] != myproc)
      ++sync_recvcounts_[proc[jj]];
    map_glb_to_loc[glb_idx[jj]] = jj;
  }

  ierr = MPI_Alltoall( &sync_recvcounts_[0], 1, MPI_INT, &sync_sendcounts_[0], 1, MPI_INT, MPI_COMM_WORLD );
  sync_recvdispls_[0]=0;
  sync_senddispls_[0]=0;
  for (int jproc=1; jproc<nproc; ++jproc)
  {
    sync_recvdispls_[jproc]=sync_recvcounts_[jproc-1]+sync_recvdispls_[jproc-1];
    sync_senddispls_[jproc]=sync_sendcounts_[jproc-1]+sync_senddispls_[jproc-1];
  }

  sync_recvcnt_ = std::accumulate(sync_recvcounts_.begin(),sync_recvcounts_.end(),0);
  sync_sendcnt_ = std::accumulate(sync_sendcounts_.begin(),sync_sendcounts_.end(),0);

  std::cout << myproc << " :  recv_cnt = " << sync_recvcnt_ << std::endl;
  std::cout << myproc << " :  send_cnt = " << sync_sendcnt_ << std::endl;
  sync_sendmap_.resize(sync_sendcnt_);
  sync_recvmap_.resize(sync_recvcnt_);

  std::vector<int> send_requests(sync_recvcnt_);
  std::vector<int> recv_requests(sync_sendcnt_);

  // Pack
  std::vector<int> cnt(nproc,0);
  for (int jj=0; jj<nb_nodes; ++jj)
  {
    if (proc[jj] != myproc)
    {
      send_requests[ sync_recvdispls_[proc[jj]]+cnt[proc[jj]] ] = glb_idx[jj];
      sync_recvmap_[ sync_recvdispls_[proc[jj]]+cnt[proc[jj]] ] = jj; 
      ++cnt[proc[jj]];
    }
  }

  ierr = MPI_Alltoallv( &send_requests[0], &sync_recvcounts_[0], &sync_recvdispls_[0], MPI_INT,
                        &recv_requests[0], &sync_sendcounts_[0], &sync_senddispls_[0], MPI_INT,
                        MPI_COMM_WORLD );

  for( int jj=1; jj<sync_sendcnt_; ++jj )
    sync_sendmap_[jj] = map_glb_to_loc[ recv_requests[jj] ];
}



template<>
inline void HaloExchange::create_mappings_impl<2,1>( 
    std::vector<int>& send_map, 
    std::vector<int>& recv_map,
    int nb_vars) const
{
  int nb_nodes = bounds_[0];
  int send_idx(0);
  int recv_idx(0);
  for (int n=0; n<sync_sendcnt_; ++n)
  {
    const int jnode = sync_sendmap_[n];
    for (int var=0; var<nb_vars; ++var)
    {
      const int varidx = var*nb_nodes;
      const int nodevaridx = jnode + varidx;
      send_map[send_idx++] = nodevaridx;
    }
  }

  for (int n=0; n<sync_recvcnt_; ++n)
  {
    const int jnode = sync_recvmap_[n];
    for (int var=0; var<nb_vars; ++var)
    {
      const int varidx = var*nb_nodes;
      const int nodevaridx = jnode + varidx;
      recv_map[recv_idx++] = nodevaridx;
    }
  }
}

/// create_mappings_impl<3,2>
template<>
inline void HaloExchange::create_mappings_impl<3,2>( 
    std::vector<int>& send_map, 
    std::vector<int>& recv_map,
    int nb_vars) const
{
  int nb_levs = bounds_[0];
  int nb_nodes = bounds_[1];
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
        case 1:
          create_mappings_impl<2,1>(send_map,recv_map,nb_vars);
          break;
        default:
          throw std::runtime_error("Not implemented");
      }
      break;
    }
    case 3:
    {
      switch (par_bound_)
      {
        case 2:
          create_mappings_impl<3,2>(send_map,recv_map,nb_vars);
          break;
        default:
          throw std::runtime_error("Not implemented");
      }
      break;
    }
    default:
      throw std::runtime_error("Not implemented");
  }
}

/////////////////////
}