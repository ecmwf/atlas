/*
 * (C) Copyright 1996-2014 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */



#ifndef Gather_hpp
#define Gather_hpp

#include <vector>
#include <stdexcept>

#include "atlas/mpl/MPL.hpp"
#include "atlas/util/Debug.hpp"
#include "atlas/util/ArrayView.hpp"

namespace atlas {

class Gather
{
public:
  Gather();
  virtual ~Gather() {}

public: // methods

  /// @brief Setup
  /// @param [in] part         List of partitions
  /// @param [in] remote_idx   List of local indices on remote partitions
  /// @param [in] base         values of remote_idx start at "base"
  /// @param [in] glb_idx      List of global indices
  /// @param [in] max_glb_idx  maximum glb index we want to gather.
  ///                          To gather everything, set to val > max value in domain
  /// @param [in] parsize      size of given lists
  void setup( const int part[],
              const int remote_idx[], const int base,
              const int glb_idx[], const int max_glb_idx,
              const int parsize );

  template <typename DATA_TYPE>
  void execute( const DATA_TYPE lfield[],
                const int lvar_strides[],
                const int lvar_extents[],
                const int lvar_rank,
                DATA_TYPE gfield[],
                const int gvar_strides[],
                const int gvar_extents[],
                const int gvar_rank ) const;

  template <typename DATA_TYPE>
  void execute( DATA_TYPE lfield[],
                DATA_TYPE gfield[],
                const int nb_vars ) const;

  template <typename DATA_TYPE, int LRANK, int GRANK>
  void execute( const ArrayView<DATA_TYPE,LRANK>& lfield,
                ArrayView<DATA_TYPE,GRANK>& gfield ) const;


  int glb_dof() const { return recvcnt_; }

private: // methods
  template< typename DATA_TYPE>
  void pack_send_buffer( const DATA_TYPE field[],
                         const int var_strides[],
                         const int var_extents[],
                         int var_rank,
                         DATA_TYPE send_buffer[] ) const;

  template< typename DATA_TYPE>
  void unpack_recv_buffer( const DATA_TYPE recv_buffer[],
                           DATA_TYPE field[],
                           const int var_strides[],
                           const int var_extents[],
                           int var_rank ) const;

  template<typename DATA_TYPE, int RANK>
  void var_info( const ArrayView<DATA_TYPE,RANK>& arr,
                 std::vector<int>& varstrides,
                 std::vector<int>& varextents ) const;

private: // data
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
  int root;

  bool is_setup_;

  int parsize_;
};

template<typename DATA_TYPE>
void Gather::execute( const DATA_TYPE lfield[],
                        const int lvar_strides[],
                        const int lvar_extents[],
                        const int lvar_rank,
                      DATA_TYPE gfield[],
                        const int gvar_strides[],
                        const int gvar_extents[],
                        const int gvar_rank ) const
{
  if( ! is_setup_ )
  {
    throw eckit::SeriousBug("Gather was not setup",Here());
  }

  int tag=1;
  int ibuf;
  int lvar_size = std::accumulate(lvar_extents,lvar_extents+lvar_rank,1,std::multiplies<int>());
  int gvar_size = std::accumulate(gvar_extents,gvar_extents+gvar_rank,1,std::multiplies<int>());
  int send_size = sendcnt_ * lvar_size;
  int recv_size = recvcnt_ * gvar_size;

  std::vector<DATA_TYPE> send_buffer(send_size);
  std::vector<DATA_TYPE> recv_buffer(recv_size);
  std::vector<MPI_Request> send_req(nproc);
  std::vector<MPI_Request> recv_req(nproc);
  std::vector<int> recv_displs(nproc);
  std::vector<int> recv_counts(nproc);


  for (int jproc=0; jproc<nproc; ++jproc)
  {
    recv_counts[jproc] = recvcounts_[jproc]*gvar_size;
    recv_displs[jproc] = recvdispls_[jproc]*gvar_size;
  }

  /// -----------------------------------------------------------
  /// With mappings and everything in place, we can now call MPI

  pack_send_buffer(lfield,lvar_strides,lvar_extents,lvar_rank,send_buffer.data());

  /// Gather
  MPL_CHECK_RESULT( MPI_Gatherv( send_buffer.data(), send_size, MPL::TYPE<DATA_TYPE>(),
                      recv_buffer.data(), const_cast<int*>(recv_counts.data()),
                      const_cast<int*>(recv_displs.data()),  MPL::TYPE<DATA_TYPE>(),
                      root, MPI_COMM_WORLD ) );

  /// Unpack
  unpack_recv_buffer(recv_buffer.data(),gfield,gvar_strides,gvar_extents,gvar_rank);
}

template<typename DATA_TYPE>
void Gather::pack_send_buffer( const DATA_TYPE field[],
                               const int var_strides[],
                               const int var_extents[],
                               int var_rank,
                               DATA_TYPE send_buffer[] ) const
{
  int ibuf = 0;
  int send_stride = var_strides[0]*var_extents[0];

  switch( var_rank )
  {
  case 1:
    for( int p=0; p<sendcnt_; ++p)
    {
      const int pp = send_stride*sendmap_[p];
      for( int i=0; i<var_extents[0]; ++i )
        send_buffer[ibuf++] = field[pp+i*var_strides[0]];
    }
    break;
  case 2:
    for( int p=0; p<sendcnt_; ++p)
    {
      const int pp = send_stride*sendmap_[p];
      for( int i=0; i<var_extents[0]; ++i )
      {
        for( int j=0; j<var_extents[1]; ++j )
        {
          send_buffer[ibuf++] = field[pp+i*var_strides[0]+j*var_strides[1]];
        }
      }
    }
    break;
  case 3:
    for( int p=0; p<sendcnt_; ++p)
    {
      const int pp = send_stride*sendmap_[p];
      for( int i=0; i<var_extents[0]; ++i )
      {
        for( int j=0; j<var_extents[1]; ++j )
        {
          for( int k=0; k<var_extents[2]; ++k )
          {
            send_buffer[ibuf++] =
              field[ pp+i*var_strides[0]+j*var_strides[1]+k*var_strides[2]];
          }
        }
      }
    }
    break;
  default:
    NOTIMP;
  }
}

template<typename DATA_TYPE>
void Gather::unpack_recv_buffer( const DATA_TYPE recv_buffer[],
                                       DATA_TYPE field[],
                                       const int var_strides[],
                                       const int var_extents[],
                                       int var_rank ) const
{
  bool field_changed = false;
  DATA_TYPE tmp;
  int ibuf = 0;
  int recv_stride = var_strides[0]*var_extents[0];

  switch( var_rank )
  {
  case 1:
    for( int p=0; p<recvcnt_; ++p)
    {
      const int pp = recv_stride*recvmap_[p];
      for( int i=0; i<var_extents[0]; ++i)
      {
        field[ pp + i*var_strides[0] ] = recv_buffer[ibuf++];
      }
    }
    break;
  case 2:
    for( int p=0; p<recvcnt_; ++p)
    {
      const int pp = recv_stride*recvmap_[p];
      for( int i=0; i<var_extents[0]; ++i )
      {
        for( int j=0; j<var_extents[1]; ++j )
        {
          field[ pp + i*var_strides[0] + j*var_strides[1] ]
              = recv_buffer[ibuf++];
        }
      }
    }
    break;
  case 3:
    for( int p=0; p<recvcnt_; ++p)
    {
      const int pp = recv_stride*recvmap_[p];
      for( int i=0; i<var_extents[0]; ++i )
      {
        for( int j=0; j<var_extents[1]; ++j )
        {
          for( int k=0; k<var_extents[2]; ++k )
          {
            field[ pp + i*var_strides[0] + j*var_strides[1] + k*var_strides[2] ]
                = recv_buffer[ibuf++];
          }
        }
      }
    }
    break;
  default:
    NOTIMP;
  }
}


template <typename DATA_TYPE>
void Gather::execute( DATA_TYPE lfield[],
                      DATA_TYPE gfield[],
                      const int nb_vars ) const
{
  int strides[] = {1};
  int extents[] = {nb_vars};
  execute( lfield, strides, extents, 1,
           gfield, strides, extents, 1 );
}


template<typename DATA_TYPE, int RANK>
void Gather::var_info( const ArrayView<DATA_TYPE,RANK>& arr,
                             std::vector<int>& varstrides,
                             std::vector<int>& varextents ) const
{
  int rank = std::max(1,RANK-1) ;
  varstrides.resize(rank);
  varextents.resize(rank);
  if( RANK>1 )
  {
    varstrides.assign(arr.strides()+1,arr.strides()+RANK);
    varextents.assign(arr.extents()+1,arr.extents()+RANK);
  }
  else
  {
    varstrides[0] = arr.strides()[0];
    varextents[0] = 1;
  }
}

template <typename DATA_TYPE, int LRANK, int GRANK>
void Gather::execute( const ArrayView<DATA_TYPE,LRANK>& lfield,
                      ArrayView<DATA_TYPE,GRANK>& gfield ) const
{
  if( lfield.size() == parsize_ && gfield.size() == recvcnt_ )
  {
    std::vector<int> lvarstrides, lvarextents;
    std::vector<int> gvarstrides, gvarextents;
    var_info(lfield, lvarstrides, lvarextents);
    var_info(gfield, gvarstrides, gvarextents);
    execute( lfield.data(), lvarstrides.data(), lvarextents.data(), lvarstrides.size(),
             gfield.data(), gvarstrides.data(), gvarextents.data(), gvarstrides.size() );
  }
  else
  {
    DEBUG_VAR(lfield.size());
    DEBUG_VAR(gfield.size());
    NOTIMP; // Need to implement with parallel ranks > 1
  }
}


// ------------------------------------------------------------------
// C wrapper interfaces to C++ routines
extern "C" 
{
  Gather* atlas__Gather__new ();
  void atlas__Gather__delete (Gather* This);
  void atlas__Gather__setup (Gather* This, int part[], int remote_idx[], int base, int glb_idx[], int max_glb_idx, int parsize);
  int atlas__Gather__glb_dof (Gather* This);
  void atlas__Gather__execute_strided_int (Gather* This, int lfield[], int lvar_strides[], int lvar_extents[], int lvar_rank, int gfield[], int gvar_strides[], int gvar_extents[], int gvar_rank);
  void atlas__Gather__execute_strided_float (Gather* This, float lfield[], int lvar_strides[], int lvar_extents[], int lvar_rank, float gfield[], int gvar_strides[], int gvar_extents[], int gvar_rank);
  void atlas__Gather__execute_strided_double (Gather* This, double lfield[], int lvar_strides[], int lvar_extents[], int lvar_rank, double gfield[], int gvar_strides[], int gvar_extents[], int gvar_rank);
  void atlas__Gather__execute_int (Gather* This, int locfield[], int glbfield[], int nb_vars);
  void atlas__Gather__execute_float (Gather* This, float locfield[], float glbfield[], int nb_vars);
  void atlas__Gather__execute_double (Gather* This, double locfield[], double glbfield[], int nb_vars);
}
// ------------------------------------------------------------------


} // namespace atlas

#endif // Gather_hpp
