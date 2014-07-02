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

class GatherScatter
{
public:
  GatherScatter();
  virtual ~GatherScatter() {}

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
  void gather( const DATA_TYPE lfield[],
               const int lvar_strides[],
               const int lvar_extents[],
               const int lvar_rank,
               DATA_TYPE gfield[],
               const int gvar_strides[],
               const int gvar_extents[],
               const int gvar_rank ) const;

  template <typename DATA_TYPE>
  void gather( const DATA_TYPE lfield[],
               DATA_TYPE gfield[],
               const int nb_vars ) const;

  template <typename DATA_TYPE, int LRANK, int GRANK>
  void gather( const ArrayView<DATA_TYPE,LRANK>& lfield,
               ArrayView<DATA_TYPE,GRANK>& gfield ) const;

  template <typename DATA_TYPE>
  void scatter( const DATA_TYPE gfield[],
                const int gvar_strides[],
                const int gvar_extents[],
                const int gvar_rank,
                DATA_TYPE lfield[],
                const int lvar_strides[],
                const int lvar_extents[],
                const int lvar_rank ) const;

  template <typename DATA_TYPE, int GRANK, int LRANK>
  void scatter( const ArrayView<DATA_TYPE,GRANK>& gfield,
                ArrayView<DATA_TYPE,LRANK>& lfield ) const;

  int glb_dof() const { return glbcnt_; }

  int loc_dof() const { return loccnt_; }

private: // methods
  template< typename DATA_TYPE>
  void pack_send_buffer( const DATA_TYPE field[],
                         const int var_strides[],
                         const int var_extents[],
                         const int var_rank,
                         const std::vector<int>& sendmap,
                         DATA_TYPE send_buffer[] ) const;

  template< typename DATA_TYPE>
  void unpack_recv_buffer( const std::vector<int>& recvmap,
                           const DATA_TYPE recv_buffer[],
                           DATA_TYPE field[],
                           const int var_strides[],
                           const int var_extents[],
                           const int var_rank ) const;

  template<typename DATA_TYPE, int RANK>
  void var_info( const ArrayView<DATA_TYPE,RANK>& arr,
                 std::vector<int>& varstrides,
                 std::vector<int>& varextents ) const;

private: // data
  int               loccnt_;
  int               glbcnt_;
  std::vector<int>  loccounts_;
  std::vector<int>  locdispls_;
  std::vector<int>  glbcounts_;
  std::vector<int>  glbdispls_;
  std::vector<int>  locmap_;
  std::vector<int>  glbmap_;

  int nproc;
  int myproc;
  int root;

  bool is_setup_;

  int parsize_;
};


////////////////////////////////////////////////////////////////////////////////


template<typename DATA_TYPE>
void GatherScatter::gather( const DATA_TYPE lfield[],
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

  const int lvar_size = std::accumulate(lvar_extents,lvar_extents+lvar_rank,1,std::multiplies<int>());
  const int gvar_size = std::accumulate(gvar_extents,gvar_extents+gvar_rank,1,std::multiplies<int>());
  const int loc_size = loccnt_ * lvar_size;
  const int glb_size = glbcnt_ * gvar_size;
  std::vector<DATA_TYPE> loc_buffer(loc_size);
  std::vector<DATA_TYPE> glb_buffer(glb_size);
  std::vector<int> glb_displs(nproc);
  std::vector<int> glb_counts(nproc);

  for (int jproc=0; jproc<nproc; ++jproc)
  {
    glb_counts[jproc] = glbcounts_[jproc]*gvar_size;
    glb_displs[jproc] = glbdispls_[jproc]*gvar_size;
  }

  /// Pack
  pack_send_buffer(lfield,lvar_strides,lvar_extents,lvar_rank,locmap_,loc_buffer.data());

  /// Gather
  MPL_CHECK_RESULT(
      MPI_Gatherv( loc_buffer.data(), loc_size, MPL::TYPE<DATA_TYPE>(),
                   glb_buffer.data(), glb_counts.data(), glb_displs.data(), MPL::TYPE<DATA_TYPE>(),
                   root, MPI_COMM_WORLD ) );

  /// Unpack
  unpack_recv_buffer(glbmap_,glb_buffer.data(),gfield,gvar_strides,gvar_extents,gvar_rank);
}



template<typename DATA_TYPE>
void GatherScatter::scatter( const DATA_TYPE gfield[],
                             const int gvar_strides[],
                             const int gvar_extents[],
                             const int gvar_rank,
                             DATA_TYPE lfield[],
                             const int lvar_strides[],
                             const int lvar_extents[],
                             const int lvar_rank ) const
{
  if( ! is_setup_ )
  {
    throw eckit::SeriousBug("GatherScatter was not setup",Here());
  }

  const int lvar_size = std::accumulate(lvar_extents,lvar_extents+lvar_rank,1,std::multiplies<int>());
  const int gvar_size = std::accumulate(gvar_extents,gvar_extents+gvar_rank,1,std::multiplies<int>());
  const int glb_size = glbcnt_ * gvar_size;
  const int loc_size = loccnt_ * lvar_size;
  std::vector<DATA_TYPE> glb_buffer(glb_size);
  std::vector<DATA_TYPE> loc_buffer(loc_size);
  std::vector<int> glb_displs(nproc);
  std::vector<int> glb_counts(nproc);

  for (int jproc=0; jproc<nproc; ++jproc)
  {
    glb_counts[jproc] = glbcounts_[jproc]*gvar_size;
    glb_displs[jproc] = glbdispls_[jproc]*gvar_size;
  }

  /// Pack
  pack_send_buffer(gfield,gvar_strides,gvar_extents,gvar_rank,glbmap_,glb_buffer.data());

  /// Scatter
  MPL_CHECK_RESULT(
      MPI_Scatterv( glb_buffer.data(), glb_counts.data(), glb_displs.data(), MPL::TYPE<DATA_TYPE>(),
                    loc_buffer.data(), loc_size, MPL::TYPE<DATA_TYPE>(),
                    root, MPI_COMM_WORLD ) );

  /// Unpack
  unpack_recv_buffer(locmap_,loc_buffer.data(),lfield,lvar_strides,lvar_extents,lvar_rank);
}

template<typename DATA_TYPE>
void GatherScatter::pack_send_buffer( const DATA_TYPE field[],
                                      const int var_strides[],
                                      const int var_extents[],
                                      const int var_rank,
                                      const std::vector<int>& sendmap,
                                      DATA_TYPE send_buffer[] ) const
{
  const int sendcnt = sendmap.size();
  int ibuf = 0;
  const int send_stride = var_strides[0]*var_extents[0];

  //DEBUG_VAR(var_rank);
  switch( var_rank )
  {
  case 1:
    for( int p=0; p<sendcnt; ++p)
    {
      const int pp = send_stride*sendmap[p];
      //DEBUG("p " << locmap_[p] << "  " << var_extents[0] << "  " << var_strides[0] << "    pp " << pp << "  ibuf " << ibuf);
      for( int i=0; i<var_extents[0]; ++i )
      {
        DATA_TYPE tmp =  field[pp+i*var_strides[0]];
        send_buffer[ibuf++] = tmp;
      }
    }
    break;
  case 2:
    for( int p=0; p<sendcnt; ++p)
    {
      const int pp = send_stride*sendmap[p];
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
    for( int p=0; p<sendcnt; ++p)
    {
      const int pp = send_stride*sendmap[p];
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
void GatherScatter::unpack_recv_buffer( const std::vector<int>& recvmap,
                                        const DATA_TYPE recv_buffer[],
                                        DATA_TYPE field[],
                                        const int var_strides[],
                                        const int var_extents[],
                                        const int var_rank ) const
{
  const int recvcnt = recvmap.size();
  bool field_changed = false;
  DATA_TYPE tmp;
  int ibuf = 0;
  const int recv_stride = var_strides[0]*var_extents[0];

  switch( var_rank )
  {
  case 1:
    for( int p=0; p<recvcnt; ++p)
    {
      const int pp = recv_stride*recvmap[p];
      for( int i=0; i<var_extents[0]; ++i)
      {
        field[ pp + i*var_strides[0] ] = recv_buffer[ibuf++];
      }
    }
    break;
  case 2:
    for( int p=0; p<recvcnt; ++p)
    {
      const int pp = recv_stride*recvmap[p];
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
    for( int p=0; p<recvcnt; ++p)
    {
      const int pp = recv_stride*recvmap[p];
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
void GatherScatter::gather( const DATA_TYPE lfield[],
                            DATA_TYPE gfield[],
                            const int nb_vars ) const
{
  int strides[] = {1};
  int extents[] = {nb_vars};
  gather( lfield, strides, extents, 1,
           gfield, strides, extents, 1 );
}


template<typename DATA_TYPE, int RANK>
void GatherScatter::var_info( const ArrayView<DATA_TYPE,RANK>& arr,
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
void GatherScatter::gather( const ArrayView<DATA_TYPE,LRANK>& lfield,
                            ArrayView<DATA_TYPE,GRANK>& gfield ) const
{
  if( lfield.size() == parsize_ && gfield.size() == glbcnt_ )
  {
    std::vector<int> lvarstrides, lvarextents;
    std::vector<int> gvarstrides, gvarextents;
    var_info(lfield, lvarstrides, lvarextents);
    var_info(gfield, gvarstrides, gvarextents);
    gather( lfield.data(), lvarstrides.data(), lvarextents.data(), lvarstrides.size(),
            gfield.data(), gvarstrides.data(), gvarextents.data(), gvarstrides.size() );
  }
  else
  {
    DEBUG_VAR(lfield.size());
    DEBUG_VAR(gfield.size());
    NOTIMP; // Need to implement with parallel ranks > 1
  }
}

template <typename DATA_TYPE, int GRANK, int LRANK>
void GatherScatter::scatter( const ArrayView<DATA_TYPE,GRANK>& gfield,
                             ArrayView<DATA_TYPE,LRANK>& lfield ) const
{
  if( lfield.size() == parsize_ && gfield.size() == glbcnt_ )
  {
    std::vector<int> lvarstrides, lvarextents;
    std::vector<int> gvarstrides, gvarextents;
    var_info(lfield, lvarstrides, lvarextents);
    var_info(gfield, gvarstrides, gvarextents);
    scatter( gfield.data(), gvarstrides.data(), gvarextents.data(), gvarstrides.size(),
             lfield.data(), lvarstrides.data(), lvarextents.data(), lvarstrides.size() );
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
  GatherScatter* atlas__GatherScatter__new ();
  void atlas__GatherScatter__delete (GatherScatter* This);
  void atlas__GatherScatter__setup (GatherScatter* This, int part[], int remote_idx[], int base, int glb_idx[], int max_glb_idx, int parsize);
  int atlas__GatherScatter__glb_dof (GatherScatter* This);
  void atlas__GatherScatter__gather_int (GatherScatter* This, int lfield[], int lvar_strides[], int lvar_extents[], int lvar_rank, int gfield[], int gvar_strides[], int gvar_extents[], int gvar_rank);
  void atlas__GatherScatter__gather_float (GatherScatter* This, float lfield[], int lvar_strides[], int lvar_extents[], int lvar_rank, float gfield[], int gvar_strides[], int gvar_extents[], int gvar_rank);
  void atlas__GatherScatter__gather_double (GatherScatter* This, double lfield[], int lvar_strides[], int lvar_extents[], int lvar_rank, double gfield[], int gvar_strides[], int gvar_extents[], int gvar_rank);
  void atlas__GatherScatter__scatter_int (GatherScatter* This, int gfield[], int gvar_strides[], int gvar_extents[], int gvar_rank, int lfield[], int lvar_strides[], int lvar_extents[], int lvar_rank);
  void atlas__GatherScatter__scatter_float (GatherScatter* This, float gfield[], int gvar_strides[], int gvar_extents[], int gvar_rank, float lfield[], int lvar_strides[], int lvar_extents[], int lvar_rank);
  void atlas__GatherScatter__scatter_double (GatherScatter* This, double gfield[], int gvar_strides[], int gvar_extents[], int gvar_rank, double lfield[], int lvar_strides[], int lvar_extents[], int lvar_rank);
}

// ------------------------------------------------------------------

//typedef GatherScatter Gather;

} // namespace atlas

#endif // Gather_hpp
