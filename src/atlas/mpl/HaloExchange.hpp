/*
 * (C) Copyright 1996-2014 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */



#ifndef HaloExchange_hpp
#define HaloExchange_hpp

#include <vector>
#include <stdexcept>

#include "eckit/exception/Exceptions.h"
#include "atlas/mpl/MPL.hpp"

namespace atlas {

class HaloExchange
{
public:
  HaloExchange();
  virtual ~HaloExchange() {}

public: // methods

  void setup( const int proc[], 
              const int glb_idx[], 
              const int periodicicity[],
              const std::vector<int>& bounds, 
              int par_bound );

  template <typename DATA_TYPE>
  void execute( DATA_TYPE field[], int nb_vars ) const;

private: // methods

  void create_mappings( std::vector<int>& send_map, std::vector<int>& recv_map, int nb_vars ) const;

  template<int N, int P>
  void create_mappings_impl( std::vector<int>& send_map, std::vector<int>& recv_map, int nb_vars ) const;

  int index(int i, int j, int k, int ni, int nj, int nk) const { return( i + ni*( j + nj*k) ); }

  int index(int i, int j, int ni, int nj) const { return( i + ni*j ); }


private: // data
  int               packet_size_;
  int               sendcnt_;
  int               recvcnt_;
  std::vector<int>  sendcounts_;
  std::vector<int>  senddispls_;
  std::vector<int>  recvcounts_;
  std::vector<int>  recvdispls_;
  std::vector<int>  sendmap_;
  std::vector<int>  recvmap_;

  int nproc;
  int myproc;

  std::vector<int> bounds_;
  int par_bound_;
  bool is_setup_;
};

template<typename DATA_TYPE>
void HaloExchange::execute( DATA_TYPE field[], int nb_vars ) const
{
  if( ! is_setup_ )
  {
    throw eckit::SeriousBug("HaloExchange was not setup",Here());
  }
#define FIELD_CONTIGUOUS true

  int tag=1;
  int ibuf;
  int point_size = packet_size_ * nb_vars;
  int send_size = sendcnt_ * point_size;
  int recv_size = recvcnt_ * point_size;

#ifndef STACK_ARRAYS
  std::vector<DATA_TYPE> send_buffer(send_size);
  std::vector<DATA_TYPE> recv_buffer(recv_size);
  std::vector<MPI_Request> send_req(nproc);
  std::vector<MPI_Request> recv_req(nproc);
  std::vector<int> send_displs(nproc);
  std::vector<int> recv_displs(nproc);
  std::vector<int> send_counts(nproc);
  std::vector<int> recv_counts(nproc);
#else // This seems to be slower on Intel ICC 13.0.1
  DATA_TYPE send_buffer[send_size];
  DATA_TYPE recv_buffer[recv_size];
  MPI_Request send_req[nproc];
  MPI_Request recv_req[nproc];
  int send_displs[nproc];
  int recv_displs[nproc];
  int send_counts[nproc];
  int recv_counts[nproc];
#endif


  // std::cout << myproc << "  :  field before = ";
  // for( int i=0; i< nb_vars*bounds_[par_bound_]; ++i)
  //   std::cout << field[i] << " ";
  // std::cout << std::endl;

  for (int jproc=0; jproc<nproc; ++jproc)
  {
    send_counts[jproc] = sendcounts_[jproc]*point_size;
    recv_counts[jproc] = recvcounts_[jproc]*point_size;
    send_displs[jproc] = senddispls_[jproc]*point_size;
    recv_displs[jproc] = recvdispls_[jproc]*point_size;
  }

#ifndef FIELD_CONTIGUOUS
  // Create additional mapping
  std::vector<int> send_map(send_size);
  std::vector<int> recv_map(recv_size);
  create_mappings(send_map,recv_map,nb_vars);
#endif
  // std::cout << myproc << "  :  send_map  = ";
  // for( int i=0; i< send_map.size(); ++i)
  //   std::cout << send_map[i] << " ";
  // std::cout << std::endl;

  // std::cout << myproc << "  :  recv_map  = ";
  // for( int i=0; i< recv_map.size(); ++i)
  //   std::cout << recv_map[i] << " ";
  // std::cout << std::endl;

  /// -----------------------------------------------------------
  /// With mappings and everything in place, we can now call MPI

  /// Let MPI know what we like to receive
  for( int jproc=0; jproc<nproc; ++jproc )
  {
    if(recv_counts[jproc] > 0)
    {
      MPL_CHECK_RESULT( MPI_Irecv( &recv_buffer[recv_displs[jproc]] , recv_counts[jproc],
        MPL::TYPE<DATA_TYPE>(), jproc, tag, MPI_COMM_WORLD, &recv_req[jproc] ) );
    }
  }

  /// Pack
#ifndef FIELD_CONTIGUOUS
  // Use additional mapping
  for( int jj=0; jj<send_size; ++jj)
    send_buffer[jj] = field[ send_map[jj] ];
#else
  // Use original mapping + contiguous bits
  ibuf = 0;
  for( int jj=0; jj<sendcnt_; ++jj)
  {
    const int ii = point_size*sendmap_[jj];
    for( int ip=0; ip<point_size; ++ip )
      send_buffer[ibuf++] = field[ ii + ip ];
  }
#endif

  /// Send
  for( int jproc=0; jproc<nproc; ++jproc )
  {
    if(send_counts[jproc] > 0)
    {
      MPL_CHECK_RESULT( MPI_Isend( &send_buffer[send_displs[jproc]], send_counts[jproc],
        MPL::TYPE<DATA_TYPE>(), jproc, tag, MPI_COMM_WORLD, &send_req[jproc] ) );
    }
  }

  /// Wait for receiving to finish
  for (int jproc=0; jproc<nproc; ++jproc)
  {
    if( recvcounts_[jproc] > 0)
    {
      MPL_CHECK_RESULT( MPI_Wait(&recv_req[jproc], MPI_STATUS_IGNORE ) );
    }
  }

  /// Unpack
#ifndef FIELD_CONTIGUOUS
  // Use additional mapping
  for( int jj=0; jj<recv_size; ++jj)
  {
    field[ recv_map[jj] ] = recv_buffer[jj];
  }
#else
  // Use original mapping + contiguous bits
  ibuf = 0;
  for( int jj=0; jj<recvcnt_; ++jj)
  {
    const int ii = point_size*recvmap_[jj];
    for( int ip=0; ip<point_size; ++ip)
      field[ ii + ip ] = recv_buffer[ibuf++];
  }
#endif

  /// Wait for sending to finish
  for (int jproc=0; jproc<nproc; ++jproc)
  {
    if( sendcounts_[jproc] > 0)
    {
      MPL_CHECK_RESULT( MPI_Wait(&send_req[jproc], MPI_STATUS_IGNORE ) );
    }
  }

  // std::cout << myproc << "  :  field after  = ";
  // for( int i=0; i< nb_vars*bounds_[par_bound_]; ++i)
  //   std::cout << field[i] << " ";
  // std::cout << std::endl;
}

// ------------------------------------------------------------------
// C wrapper interfaces to C++ routines
extern "C" 
{
  HaloExchange* atlas__HaloExchange__new ();
  void atlas__HaloExchange__delete (HaloExchange* This);
  void atlas__HaloExchange__setup (HaloExchange* This, int proc[], int glb_idx[], int master_glb_idx[], int bounds[], int nb_bounds, int par_bound);
  void atlas__HaloExchange__execute_int (HaloExchange* This, int field[], int nb_vars);
  void atlas__HaloExchange__execute_float (HaloExchange* This, float field[], int nb_vars);
  void atlas__HaloExchange__execute_double (HaloExchange* This, double field[], int nb_vars);

}
// ------------------------------------------------------------------


} // namespace atlas

#endif // HaloExchange_hpp
