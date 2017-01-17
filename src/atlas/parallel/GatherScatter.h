/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef Gather_h
#define Gather_h

#include <vector>
#include <stdexcept>

#include "eckit/memory/SharedPtr.h"
#include "eckit/memory/Owned.h"
#include "atlas/parallel/mpi/mpi.h"

#include "atlas/internals/atlas_config.h"
#include "atlas/internals/Debug.h"
#include "atlas/array/ArrayView.h"
#include "atlas/internals/MPLArrayView.h"

namespace atlas {
namespace parallel {

//----------------------------------------------------------------------------------------------------------------------

template<typename T> struct remove_const          { typedef T type; };
template<typename T> struct remove_const<T const> { typedef T type; };

template<typename T> struct add_const          { typedef const typename remove_const<T>::type type; };
template<typename T> struct add_const<T const> { typedef const T type; };

template <typename DATA_TYPE>
class Field
{
private:
  typedef typename remove_const<DATA_TYPE>::type     NON_CONST_DATA_TYPE;
  typedef typename add_const<DATA_TYPE>::type        CONST_DATA_TYPE;

public:
  Field() {}

  Field( DATA_TYPE data_[], const size_t var_strides_[], const size_t var_shape_[], const size_t var_rank_ ) :
    data(const_cast<DATA_TYPE*>(data_)), var_rank(var_rank_)
  {
    var_strides.assign(var_strides_,var_strides_+var_rank_);
    var_shape.assign(var_shape_,var_shape_+var_rank_);
  }

  Field( DATA_TYPE data_[], const size_t nb_vars ) :
    data(const_cast<DATA_TYPE*>(data_)), var_rank(1)
  {
    var_strides.assign(1,1);
    var_shape.assign(1,nb_vars);
  }

  template<int RANK>
  Field( const array::ArrayView<NON_CONST_DATA_TYPE,RANK>& arr )
  {
    data = const_cast<DATA_TYPE*>(arr.data());
    var_rank = std::max(1,(int)arr.rank()-1);
    var_strides.resize(var_rank);
    var_shape.resize(var_rank);
    if( arr.rank()>1 )
    {
      var_strides.assign(arr.strides()+1,arr.strides()+arr.rank());
      var_shape.assign(arr.shape()+1,arr.shape()+arr.rank());
    }
    if( arr.rank() == 1 )
    {
      var_strides[0] = arr.strides()[0];
      var_shape[0] = 1;
    }
  }

public:
  DATA_TYPE* data;
  std::vector<size_t> var_strides;
  std::vector<size_t> var_shape;
  size_t var_rank;
};

class GatherScatter: public eckit::Owned {

public: // types
  typedef eckit::SharedPtr<GatherScatter> Ptr;
public:
  GatherScatter();
  GatherScatter(const std::string& name);
  virtual ~GatherScatter() {}

public: // methods

  std::string name() const { return name_; }

  /// @brief Setup
  /// @param [in] part         List of partitions
  /// @param [in] remote_idx   List of local indices on remote partitions
  /// @param [in] base         values of remote_idx start at "base"
  /// @param [in] glb_idx      List of global indices
  /// @param [in] max_glb_idx  maximum glb index we want to gather.
  ///                          To gather everything, set to val > max value in domain
  /// @param [in] parsize      size of given lists
  /// @param [in] include_ghost    Warning: setting this to true breaks scatter
  ///                              functionality. This feature allows to check periodic
  ///                              halo values in output when TRUE
  void setup( const int part[],
              const int remote_idx[], const int base,
              const gidx_t glb_idx[],
              const size_t parsize );


  /// @brief Setup
  /// @param [in] part         List of partitions
  /// @param [in] remote_idx   List of local indices on remote partitions
  /// @param [in] base         values of remote_idx start at "base"
  /// @param [in] glb_idx      List of global indices
  /// @param [in] parsize      size of given lists
  /// @param [in] mask         Mask indices not to include in the communication
  ///                          pattern (0=include,1=exclude)
  void setup( const int part[],
              const int remote_idx[], const int base,
              const gidx_t glb_idx[], const int mask[],  const size_t parsize );

  template <typename DATA_TYPE>
  void gather( const DATA_TYPE ldata[],
               const size_t lstrides[],
               const size_t lshape[],
               const size_t lrank,
               const size_t lmpl_idxpos[],
               const size_t lmpl_rank,
               DATA_TYPE gdata[],
               const size_t gstrides[],
               const size_t gshape[],
               const size_t grank,
               const size_t gmpl_idxpos[],
               const size_t gmpl_rank,
               const size_t root ) const;

  template <typename DATA_TYPE>
  void gather( const DATA_TYPE ldata[],
               const size_t lvar_strides[],
               const size_t lvar_shape[],
               const size_t lvar_rank,
               DATA_TYPE gdata[],
               const size_t gvar_strides[],
               const size_t gvar_shape[],
               const size_t gvar_rank,
               const size_t root = 0 ) const;

  template <typename DATA_TYPE>
  void gather( parallel::Field<DATA_TYPE const> lfields[],
               parallel::Field<DATA_TYPE      > gfields[],
               const size_t nb_fields,
               const size_t root = 0 ) const;

  template <typename DATA_TYPE, int LRANK, int GRANK>
  void gather( const array::ArrayView<DATA_TYPE,LRANK>& ldata,
               array::ArrayView<DATA_TYPE,GRANK>& gdata,
               const size_t root = 0 ) const;

  template <typename DATA_TYPE>
  void scatter( parallel::Field<DATA_TYPE const> gfields[],
                parallel::Field<DATA_TYPE      > lfields[],
                const size_t nb_fields,
                const size_t root = 0 ) const;

  template <typename DATA_TYPE>
  void scatter( const DATA_TYPE gdata[],
                const size_t gstrides[],
                const size_t gshape[],
                const size_t grank,
                const size_t gmpl_idxpos[],
                const size_t gmpl_rank,
                DATA_TYPE ldata[],
                const size_t lstrides[],
                const size_t lshape[],
                const size_t lrank,
                const size_t lmpl_idxpos[],
                const size_t lmpl_rank,
                const size_t root ) const;

  template <typename DATA_TYPE>
  void scatter( const DATA_TYPE gdata[],
                const size_t gvar_strides[],
                const size_t gvar_shape[],
                const size_t gvar_rank,
                DATA_TYPE ldata[],
                const size_t lvar_strides[],
                const size_t lvar_shape[],
                const size_t lvar_rank,
                const size_t root = 0 ) const;

  template <typename DATA_TYPE, int GRANK, int LRANK>
  void scatter( const array::ArrayView<DATA_TYPE,GRANK>& gdata,
                array::ArrayView<DATA_TYPE,LRANK>& ldata,
                const size_t root = 0 ) const;

  int glb_dof() const { return glbcnt_; }

  int loc_dof() const { return loccnt_; }

private: // methods
  template< typename DATA_TYPE>
  void pack_send_buffer( const parallel::Field<DATA_TYPE const>& field,
                         const std::vector<int>& sendmap,
                         DATA_TYPE send_buffer[] ) const;

  template< typename DATA_TYPE>
  void unpack_recv_buffer( const std::vector<int>& recvmap,
                           const DATA_TYPE recv_buffer[],
                           const parallel::Field<DATA_TYPE>& field ) const;

  template<typename DATA_TYPE, int RANK>
  void var_info( const array::ArrayView<DATA_TYPE,RANK>& arr,
                 std::vector<size_t>& varstrides,
                 std::vector<size_t>& varshape ) const;

private: // data

  std::string name_;
  int         loccnt_;
  int         glbcnt_;
  std::vector<int>  glbcounts_;
  std::vector<int>  glbdispls_;
  std::vector<int>  locmap_;
  std::vector<int>  glbmap_;

  size_t nproc;
  size_t myproc;

  bool is_setup_;

  size_t parsize_;


  int glb_cnt(size_t root) const { return myproc==root ? glbcnt_ : 0 ; }
};


////////////////////////////////////////////////////////////////////////////////

template <typename DATA_TYPE>
void GatherScatter::gather( parallel::Field<DATA_TYPE const> lfields[],
                            parallel::Field<DATA_TYPE> gfields[],
                            size_t nb_fields,
                            const size_t root ) const
{
  if( ! is_setup_ )
  {
    throw eckit::SeriousBug("GatherScatter was not setup",Here());
  }

  for( size_t jfield=0; jfield<nb_fields; ++jfield )
  {
    const size_t lvar_size = std::accumulate(lfields[jfield].var_shape.data(),lfields[jfield].var_shape.data()+lfields[jfield].var_rank,1,std::multiplies<size_t>());
    const size_t gvar_size = std::accumulate(gfields[jfield].var_shape.data(),gfields[jfield].var_shape.data()+gfields[jfield].var_rank,1,std::multiplies<size_t>());
    const int loc_size = loccnt_ * lvar_size;
    const int glb_size = glb_cnt(root) * gvar_size;
    std::vector<DATA_TYPE> loc_buffer(loc_size);
    std::vector<DATA_TYPE> glb_buffer(glb_size);
    std::vector<int> glb_displs(nproc);
    std::vector<int> glb_counts(nproc);

    for (size_t jproc=0; jproc<nproc; ++jproc)
    {
      glb_counts[jproc] = glbcounts_[jproc]*gvar_size;
      glb_displs[jproc] = glbdispls_[jproc]*gvar_size;
    }

    /// Pack

    pack_send_buffer(lfields[jfield],locmap_,loc_buffer.data());

    /// Gather

    parallel::mpi::comm().gatherv(loc_buffer, glb_buffer, glb_counts, glb_displs, root);

    /// Unpack
    if( myproc == root )
      unpack_recv_buffer(glbmap_,glb_buffer.data(),gfields[jfield]);
  }
}

template<typename DATA_TYPE>
void GatherScatter::gather( const DATA_TYPE ldata[],
                            const size_t lstrides[],
                            const size_t lshape[],
                            const size_t lrank,
                            const size_t lmpl_idxpos[],
                            const size_t lmpl_rank,
                            DATA_TYPE gdata[],
                            const size_t gstrides[],
                            const size_t gshape[],
                            const size_t grank,
                            const size_t gmpl_idxpos[],
                            const size_t gmpl_rank,
                            const size_t root ) const
{
  // compatibility mode
  internals::MPL_ArrayView<DATA_TYPE const> lview(ldata,lstrides,lshape,lrank,lmpl_idxpos,lmpl_rank);
  internals::MPL_ArrayView<DATA_TYPE      > gview(gdata,gstrides,gshape,grank,gmpl_idxpos,gmpl_rank);

  gather(lview.data(),lview.var_strides().data(),lview.var_shape().data(),lview.var_rank(),
         gview.data(),gview.var_strides().data(),gview.var_shape().data(),gview.var_rank(),
         root);
}

template<typename DATA_TYPE>
void GatherScatter::gather( const DATA_TYPE ldata[],
                            const size_t lvar_strides[],
                            const size_t lvar_shape[],
                            const size_t lvar_rank,
                            DATA_TYPE gdata[],
                            const size_t gvar_strides[],
                            const size_t gvar_shape[],
                            const size_t gvar_rank,
                            const size_t root ) const
{
  parallel::Field<DATA_TYPE const> lfield(ldata,lvar_strides,lvar_shape,lvar_rank);
  parallel::Field<DATA_TYPE      > gfield(gdata,gvar_strides,gvar_shape,gvar_rank);
  gather( &lfield, &gfield, 1, root );
}


template <typename DATA_TYPE>
void GatherScatter::scatter( parallel::Field<DATA_TYPE const> gfields[],
                             parallel::Field<DATA_TYPE      > lfields[],
                             const size_t nb_fields,
                             const size_t root ) const
{
  if( ! is_setup_ )
  {
    throw eckit::SeriousBug("GatherScatter was not setup",Here());
  }

  for( int jfield=0; jfield<nb_fields; ++jfield )
  {
    const int lvar_size = std::accumulate(lfields[jfield].var_shape.data(),lfields[jfield].var_shape.data()+lfields[jfield].var_rank,1,std::multiplies<int>());
    const int gvar_size = std::accumulate(gfields[jfield].var_shape.data(),gfields[jfield].var_shape.data()+gfields[jfield].var_rank,1,std::multiplies<int>());
    const int loc_size = loccnt_ * lvar_size;
    const int glb_size = glb_cnt(root) * gvar_size;
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
    if( myproc == root )
      pack_send_buffer(gfields[jfield],glbmap_,glb_buffer.data());

    /// Scatter

    parallel::mpi::comm().scatterv(glb_buffer.begin(), glb_buffer.end(), glb_counts, glb_displs, loc_buffer.begin(), loc_buffer.end(), root);

    /// Unpack
    unpack_recv_buffer(locmap_,loc_buffer.data(),lfields[jfield]);
  }
}

template<typename DATA_TYPE>
void GatherScatter::scatter( const DATA_TYPE gdata[],
                             const size_t gstrides[],
                             const size_t gshape[],
                             const size_t grank,
                             const size_t gmpl_idxpos[],
                             const size_t gmpl_rank,
                             DATA_TYPE ldata[],
                             const size_t lstrides[],
                             const size_t lshape[],
                             const size_t lrank,
                             const size_t lmpl_idxpos[],
                             const size_t lmpl_rank,
                             const size_t root ) const
{
  // compatibility mode
  internals::MPL_ArrayView<DATA_TYPE const> gview(gdata,gstrides,gshape,grank,gmpl_idxpos,gmpl_rank);
  internals::MPL_ArrayView<DATA_TYPE      > lview(ldata,lstrides,lshape,lrank,lmpl_idxpos,lmpl_rank);
  std::vector<size_t> gvar_strides(gmpl_rank);
  std::vector<size_t> gvar_shape(gmpl_rank);
  std::vector<size_t> lvar_strides(lmpl_rank);
  std::vector<size_t> lvar_shape(lmpl_rank);
  for( size_t i=0; i<gmpl_rank; ++i)
  {
    gvar_strides[i] = gview.stride(i);
    gvar_shape[i] = gview.extent(i);
  }
  for( size_t i=0; i<lmpl_rank; ++i)
  {
    lvar_strides[i] = lview.stride(i);
    lvar_shape[i] = lview.extent(i);
  }
  parallel::Field<DATA_TYPE const> gfield(gdata,gvar_strides.data(),gvar_shape.data(),gview.var_rank());
  parallel::Field<DATA_TYPE      > lfield(ldata,lvar_strides.data(),lvar_shape.data(),lview.var_rank());
  scatter( &gfield, &lfield, 1, root );
}

template<typename DATA_TYPE>
void GatherScatter::scatter( const DATA_TYPE gdata[],
                             const size_t gvar_strides[],
                             const size_t gvar_shape[],
                             const size_t gvar_rank,
                             DATA_TYPE ldata[],
                             const size_t lvar_strides[],
                             const size_t lvar_shape[],
                             const size_t lvar_rank,
                             const size_t root ) const
{
  parallel::Field<DATA_TYPE const> gfield(gdata,gvar_strides,gvar_shape,gvar_rank);
  parallel::Field<DATA_TYPE      > lfield(ldata,lvar_strides,lvar_shape,lvar_rank);
  scatter( &gfield, &lfield, 1, root );
}

template<typename DATA_TYPE>
void GatherScatter::pack_send_buffer( const parallel::Field<DATA_TYPE const>& field,
                                      const std::vector<int>& sendmap,
                                      DATA_TYPE send_buffer[] ) const
{
  const int sendcnt = sendmap.size();

  size_t ibuf = 0;
  const size_t send_stride = field.var_strides[0]*field.var_shape[0];

  switch( field.var_rank )
  {
  case 1:
    for( size_t p=0; p<sendcnt; ++p)
    {
      const size_t pp = send_stride*sendmap[p];
      for( size_t i=0; i<field.var_shape[0]; ++i )
      {
        DATA_TYPE tmp =  field.data[pp+i*field.var_strides[0]];
        send_buffer[ibuf++] = tmp;
      }
    }
    break;
  case 2:
    for( size_t p=0; p<sendcnt; ++p)
    {
      const size_t pp = send_stride*sendmap[p];
      for( size_t i=0; i<field.var_shape[0]; ++i )
      {
        const size_t ii = pp + i*field.var_strides[0];
        for( size_t j=0; j<field.var_shape[1]; ++j )
        {
          send_buffer[ibuf++] = field.data[ii+j*field.var_strides[1]];
        }
      }
    }
    break;
  case 3:
    for( size_t p=0; p<sendcnt; ++p)
    {
      const size_t pp = send_stride*sendmap[p];
      for( size_t i=0; i<field.var_shape[0]; ++i )
      {
        const size_t ii = pp + i*field.var_strides[0];
        for( size_t j=0; j<field.var_shape[1]; ++j )
        {
          const size_t jj = ii + j*field.var_strides[1];
          for( size_t k=0; k<field.var_shape[2]; ++k )
          {
            send_buffer[ibuf++] = field.data[ jj+k*field.var_strides[2]];
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
                                        const parallel::Field<DATA_TYPE>& field ) const
{
  const size_t recvcnt = recvmap.size();

  int ibuf = 0;
  const size_t recv_stride = field.var_strides[0]*field.var_shape[0];

  switch( field.var_rank )
  {
  case 1:
    for( size_t p=0; p<recvcnt; ++p)
    {
      const size_t pp = recv_stride*recvmap[p];
      for( size_t i=0; i<field.var_shape[0]; ++i)
      {
        field.data[ pp + i*field.var_strides[0] ] = recv_buffer[ibuf++];
      }
    }
    break;
  case 2:
    for( size_t p=0; p<recvcnt; ++p)
    {
      const size_t pp = recv_stride*recvmap[p];
      for( size_t i=0; i<field.var_shape[0]; ++i )
      {
        const size_t ii = pp + i*field.var_strides[0];
        for( size_t j=0; j<field.var_shape[1]; ++j )
        {
          field.data[ ii + j*field.var_strides[1] ] = recv_buffer[ibuf++];
        }
      }
    }
    break;
  case 3:
    for( size_t p=0; p<recvcnt; ++p)
    {
      const size_t pp = recv_stride*recvmap[p];
      for( size_t i=0; i<field.var_shape[0]; ++i )
      {
        const size_t ii = pp + i*field.var_strides[0];
        for( size_t j=0; j<field.var_shape[1]; ++j )
        {
          const size_t jj = ii + j*field.var_strides[1];
          for( size_t k=0; k<field.var_shape[2]; ++k )
          {
            field.data[ jj + k*field.var_strides[2] ] = recv_buffer[ibuf++];
          }
        }
      }
    }
    break;
  default:
    NOTIMP;
  }
}



template<typename DATA_TYPE, int RANK>
void GatherScatter::var_info( const array::ArrayView<DATA_TYPE,RANK>& arr,
                              std::vector<size_t>& varstrides,
                              std::vector<size_t>& varshape ) const
{
  size_t rank = std::max(1,RANK-1) ;
  varstrides.resize(rank);
  varshape.resize(rank);
  if( RANK>1 )
  {
    varstrides.assign(arr.strides()+1,arr.strides()+RANK);
    varshape.assign(arr.shape()+1,arr.shape()+RANK);
  }
  else
  {
    varstrides[0] = arr.strides()[0];
    varshape[0] = 1;
  }
}

template <typename DATA_TYPE, int LRANK, int GRANK>
void GatherScatter::gather( const array::ArrayView<DATA_TYPE,LRANK>& ldata,
                            array::ArrayView<DATA_TYPE,GRANK>& gdata,
                            const size_t root ) const
{
  if( ldata.shape(0) == parsize_ && gdata.shape(0) == glb_cnt(root) )
  {
    std::vector< parallel::Field<DATA_TYPE const> > lfields(1, parallel::Field<DATA_TYPE const>(ldata) );
    std::vector< parallel::Field<DATA_TYPE> >       gfields(1, parallel::Field<DATA_TYPE>(gdata) );
    gather( lfields.data(), gfields.data(), 1, root );
  }
  else
  {
    DEBUG_VAR(parsize_);
    DEBUG_VAR(ldata.shape(0));
    DEBUG_VAR(glb_cnt(root));
    DEBUG_VAR(gdata.shape(0));
    NOTIMP; // Need to implement with parallel ranks > 1
  }
}

template <typename DATA_TYPE, int GRANK, int LRANK>
void GatherScatter::scatter( const array::ArrayView<DATA_TYPE,GRANK>& gdata,
                             array::ArrayView<DATA_TYPE,LRANK>& ldata,
                             const size_t root ) const
{
  if( ldata.shape(0) == parsize_ && gdata.shape(0) == glb_cnt(root) )
  {
    std::vector< parallel::Field<DATA_TYPE const> > gfields(1, parallel::Field<DATA_TYPE const>(gdata) );
    std::vector< parallel::Field<DATA_TYPE> >       lfields(1, parallel::Field<DATA_TYPE>(ldata) );
    scatter( gfields.data(), lfields.data(), 1, root );
  }
  else
  {
    DEBUG_VAR(parsize_);
    DEBUG_VAR(ldata.shape(0));
    DEBUG_VAR(glb_cnt(root));
    DEBUG_VAR(gdata.shape(0));
    NOTIMP; // Need to implement with parallel ranks > 1
  }
}


// ------------------------------------------------------------------
// C wrapper interfaces to C++ routines
extern "C"
{
  GatherScatter* atlas__GatherScatter__new ();
  void atlas__GatherScatter__delete (GatherScatter* This);
  void atlas__GatherScatter__setup32 (GatherScatter* This, int part[], int remote_idx[], int base, int glb_idx[],  int parsize);
  void atlas__GatherScatter__setup64 (GatherScatter* This, int part[], int remote_idx[], int base, long glb_idx[], int parsize);
  int atlas__GatherScatter__glb_dof (GatherScatter* This);
  void atlas__GatherScatter__gather_int (GatherScatter* This, int ldata[], int lvar_strides[], int lvar_shape[], int lvar_rank, int gdata[], int gvar_strides[], int gvar_shape[], int gvar_rank);
  void atlas__GatherScatter__gather_long (GatherScatter* This, long ldata[], int lvar_strides[], int lvar_shape[], int lvar_rank, long gdata[], int gvar_strides[], int gvar_shape[], int gvar_rank);
  void atlas__GatherScatter__gather_float (GatherScatter* This, float ldata[], int lvar_strides[], int lvar_shape[], int lvar_rank, float gdata[], int gvar_strides[], int gvar_shape[], int gvar_rank);
  void atlas__GatherScatter__gather_double (GatherScatter* This, double ldata[], int lvar_strides[], int lvar_shape[], int lvar_rank, double gdata[], int gvar_strides[], int gvar_shape[], int gvar_rank);
  void atlas__GatherScatter__scatter_int (GatherScatter* This, int gdata[], int gvar_strides[], int gvar_shape[], int gvar_rank, int ldata[], int lvar_strides[], int lvar_shape[], int lvar_rank);
  void atlas__GatherScatter__scatter_long (GatherScatter* This, long gdata[], int gvar_strides[], int gvar_shape[], int gvar_rank, long ldata[], int lvar_strides[], int lvar_shape[], int lvar_rank);
  void atlas__GatherScatter__scatter_float (GatherScatter* This, float gdata[], int gvar_strides[], int gvar_shape[], int gvar_rank, float ldata[], int lvar_strides[], int lvar_shape[], int lvar_rank);
  void atlas__GatherScatter__scatter_double (GatherScatter* This, double gdata[], int gvar_strides[], int gvar_shape[], int gvar_rank, double ldata[], int lvar_strides[], int lvar_shape[], int lvar_rank);
}

// ------------------------------------------------------------------

//typedef GatherScatter Gather;

} // namespace parallel
} // namespace atlas

#endif // Gather_h
