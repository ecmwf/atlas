/*
 * (C) Copyright 1996-2014 ECMWF.
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

#include "atlas/mpl/MPL.h"
#include "atlas/util/Debug.h"
#include "atlas/util/ArrayView.h"
#include "atlas/mpl/MPLArrayView.h"

namespace atlas {

namespace mpl {
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

  Field( DATA_TYPE data_[], const int var_strides_[], const int var_extents_[], const int var_rank_ ) :
    data(const_cast<DATA_TYPE*>(data_)), var_rank(var_rank_)
  {
    var_strides.assign(var_strides_,var_strides_+var_rank_);
    var_extents.assign(var_extents_,var_extents_+var_rank_);
  }

  Field( DATA_TYPE data_[], const int nb_vars ) :
    data(const_cast<DATA_TYPE*>(data_)), var_rank(1)
  {
    var_strides.assign(1,1);
    var_extents.assign(1,nb_vars);
  }

  template<int RANK>
  Field( const ArrayView<NON_CONST_DATA_TYPE,RANK>& arr )
  {
    data = const_cast<DATA_TYPE*>(arr.data());
    var_rank = std::max(1,RANK-1);
    var_strides.resize(var_rank);
    var_extents.resize(var_rank);
    if( RANK>1 )
    {
      var_strides.assign(arr.strides()+1,arr.strides()+RANK);
      var_extents.assign(arr.extents()+1,arr.extents()+RANK);
    }
    else
    {
      var_strides[0] = arr.strides()[0];
      var_extents[0] = 1;
    }
  }

public:
  DATA_TYPE* data;
  std::vector<int> var_strides;
  std::vector<int> var_extents;
  int var_rank;
};

class GatherScatter: public eckit::Owned {

public: // types
    typedef eckit::SharedPtr<GatherScatter> Ptr;
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
  void gather( const DATA_TYPE ldata[],
               const int lstrides[],
               const int lextents[],
               const int lrank,
               const int lmpl_idxpos[],
               const int lmpl_rank,
               DATA_TYPE gdata[],
               const int gstrides[],
               const int gextents[],
               const int grank,
               const int gmpl_idxpos[],
               const int gmpl_rank,
               const int root ) const;

  template <typename DATA_TYPE>
  void gather( const DATA_TYPE ldata[],
               const int lvar_strides[],
               const int lvar_extents[],
               const int lvar_rank,
               DATA_TYPE gdata[],
               const int gvar_strides[],
               const int gvar_extents[],
               const int gvar_rank,
               const int root = 0 ) const;

  template <typename DATA_TYPE>
  void gather( mpl::Field<DATA_TYPE const> lfields[],
               mpl::Field<DATA_TYPE      > gfields[],
               const int nb_fields,
               const int root = 0 ) const;

  template <typename DATA_TYPE, int LRANK, int GRANK>
  void gather( const ArrayView<DATA_TYPE,LRANK>& ldata,
               ArrayView<DATA_TYPE,GRANK>& gdata,
               const int root = 0 ) const;

  template <typename DATA_TYPE>
  void scatter( mpl::Field<DATA_TYPE const> gfields[],
                mpl::Field<DATA_TYPE      > lfields[],
                const int nb_fields,
                const int root = 0 ) const;

  template <typename DATA_TYPE>
  void scatter( const DATA_TYPE gdata[],
                const int gstrides[],
                const int gextents[],
                const int grank,
                const int gmpl_idxpos[],
                const int gmpl_rank,
                DATA_TYPE ldata[],
                const int lstrides[],
                const int lextents[],
                const int lrank,
                const int lmpl_idxpos[],
                const int lmpl_rank,
                const int root ) const;

  template <typename DATA_TYPE>
  void scatter( const DATA_TYPE gdata[],
                const int gvar_strides[],
                const int gvar_extents[],
                const int gvar_rank,
                DATA_TYPE ldata[],
                const int lvar_strides[],
                const int lvar_extents[],
                const int lvar_rank,
                const int root = 0 ) const;

  template <typename DATA_TYPE, int GRANK, int LRANK>
  void scatter( const ArrayView<DATA_TYPE,GRANK>& gdata,
                ArrayView<DATA_TYPE,LRANK>& ldata,
                const int root = 0 ) const;

  int glb_dof() const { return glbcnt_; }

  int loc_dof() const { return loccnt_; }

private: // methods
  template< typename DATA_TYPE>
  void pack_send_buffer( const mpl::Field<DATA_TYPE const>& field,
                         const std::vector<int>& sendmap,
                         DATA_TYPE send_buffer[] ) const;

  template< typename DATA_TYPE>
  void unpack_recv_buffer( const std::vector<int>& recvmap,
                           const DATA_TYPE recv_buffer[],
                           const mpl::Field<DATA_TYPE>& field ) const;

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
  int root_;

  bool is_setup_;

  int parsize_;
};


////////////////////////////////////////////////////////////////////////////////

template <typename DATA_TYPE>
void GatherScatter::gather( mpl::Field<DATA_TYPE const> lfields[],
                            mpl::Field<DATA_TYPE> gfields[],
                            int nb_fields,
                            const int root ) const
{
  if( ! is_setup_ )
  {
    throw eckit::SeriousBug("GatherScatter was not setup",Here());
  }

  for( int jfield=0; jfield<nb_fields; ++jfield )
  {
    const int lvar_size = std::accumulate(lfields[jfield].var_extents.data(),lfields[jfield].var_extents.data()+lfields[jfield].var_rank,1,std::multiplies<int>());
    const int gvar_size = std::accumulate(gfields[jfield].var_extents.data(),gfields[jfield].var_extents.data()+gfields[jfield].var_rank,1,std::multiplies<int>());
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
    pack_send_buffer(lfields[jfield],locmap_,loc_buffer.data());

    /// Gather
    MPL_CHECK_RESULT(
        MPI_Gatherv( loc_buffer.data(), loc_size, MPL::TYPE<DATA_TYPE>(),
                     glb_buffer.data(), glb_counts.data(), glb_displs.data(), MPL::TYPE<DATA_TYPE>(),
                     root, MPI_COMM_WORLD ) );

    /// Unpack
    unpack_recv_buffer(glbmap_,glb_buffer.data(),gfields[jfield]);
  }
}

template<typename DATA_TYPE>
void GatherScatter::gather( const DATA_TYPE ldata[],
                            const int lstrides[],
                            const int lextents[],
                            const int lrank,
                            const int lmpl_idxpos[],
                            const int lmpl_rank,
                            DATA_TYPE gdata[],
                            const int gstrides[],
                            const int gextents[],
                            const int grank,
                            const int gmpl_idxpos[],
                            const int gmpl_rank,
                            const int root ) const
{
  // compatibility mode
  mpl::MPL_ArrayView<DATA_TYPE const> lview(ldata,lstrides,lextents,lrank,lmpl_idxpos,lmpl_rank);
  mpl::MPL_ArrayView<DATA_TYPE      > gview(gdata,gstrides,gextents,grank,gmpl_idxpos,gmpl_rank);

  gather(lview.data(),lview.var_strides().data(),lview.var_extents().data(),lview.var_rank(),
         gview.data(),gview.var_strides().data(),gview.var_extents().data(),gview.var_rank(),
         root);
}

template<typename DATA_TYPE>
void GatherScatter::gather( const DATA_TYPE ldata[],
                            const int lvar_strides[],
                            const int lvar_extents[],
                            const int lvar_rank,
                            DATA_TYPE gdata[],
                            const int gvar_strides[],
                            const int gvar_extents[],
                            const int gvar_rank,
                            const int root ) const
{
  mpl::Field<DATA_TYPE const> lfield(ldata,lvar_strides,lvar_extents,lvar_rank);
  mpl::Field<DATA_TYPE      > gfield(gdata,gvar_strides,gvar_extents,gvar_rank);
  gather( &lfield, &gfield, 1, root );
}


template <typename DATA_TYPE>
void GatherScatter::scatter( mpl::Field<DATA_TYPE const> gfields[],
                             mpl::Field<DATA_TYPE      > lfields[],
                             const int nb_fields,
                             const int root ) const
{
  if( ! is_setup_ )
  {
    throw eckit::SeriousBug("GatherScatter was not setup",Here());
  }

  for( int jfield=0; jfield<nb_fields; ++jfield )
  {
    const int lvar_size = std::accumulate(lfields[jfield].var_extents.data(),lfields[jfield].var_extents.data()+lfields[jfield].var_rank,1,std::multiplies<int>());
    const int gvar_size = std::accumulate(gfields[jfield].var_extents.data(),gfields[jfield].var_extents.data()+gfields[jfield].var_rank,1,std::multiplies<int>());
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
    pack_send_buffer(gfields[jfield],glbmap_,glb_buffer.data());

    /// Scatter
    MPL_CHECK_RESULT(
        MPI_Scatterv( glb_buffer.data(), glb_counts.data(), glb_displs.data(), MPL::TYPE<DATA_TYPE>(),
                      loc_buffer.data(), loc_size, MPL::TYPE<DATA_TYPE>(),
                      root, MPI_COMM_WORLD ) );

    /// Unpack
    unpack_recv_buffer(locmap_,loc_buffer.data(),lfields[jfield]);
  }
}

template<typename DATA_TYPE>
void GatherScatter::scatter( const DATA_TYPE gdata[],
                             const int gstrides[],
                             const int gextents[],
                             const int grank,
                             const int gmpl_idxpos[],
                             const int gmpl_rank,
                             DATA_TYPE ldata[],
                             const int lstrides[],
                             const int lextents[],
                             const int lrank,
                             const int lmpl_idxpos[],
                             const int lmpl_rank,
                             const int root ) const
{
  // compatibility mode
  mpl::MPL_ArrayView<DATA_TYPE const> gview(gdata,gstrides,gextents,grank,gmpl_idxpos,gmpl_rank);
  mpl::MPL_ArrayView<DATA_TYPE      > lview(ldata,lstrides,lextents,lrank,lmpl_idxpos,lmpl_rank);
  std::vector<int> gvar_strides(gmpl_rank);
  std::vector<int> gvar_extents(gmpl_rank);
  std::vector<int> lvar_strides(lmpl_rank);
  std::vector<int> lvar_extents(lmpl_rank);
  for( int i=0; i<gmpl_rank; ++i)
  {
    gvar_strides[i] = gview.stride(i);
    gvar_extents[i] = gview.extent(i);
  }
  for( int i=0; i<lmpl_rank; ++i)
  {
    lvar_strides[i] = lview.stride(i);
    lvar_extents[i] = lview.extent(i);
  }
  mpl::Field<DATA_TYPE const> gfield(gdata,gvar_strides.data(),gvar_extents.data(),gview.var_rank());
  mpl::Field<DATA_TYPE      > lfield(ldata,lvar_strides.data(),lvar_extents.data(),lview.var_rank());
  scatter( &gfield, &lfield, 1, root );
}

template<typename DATA_TYPE>
void GatherScatter::scatter( const DATA_TYPE gdata[],
                             const int gvar_strides[],
                             const int gvar_extents[],
                             const int gvar_rank,
                             DATA_TYPE ldata[],
                             const int lvar_strides[],
                             const int lvar_extents[],
                             const int lvar_rank,
                             const int root ) const
{
  mpl::Field<DATA_TYPE const> gfield(gdata,gvar_strides,gvar_extents,gvar_rank);
  mpl::Field<DATA_TYPE      > lfield(ldata,lvar_strides,lvar_extents,lvar_rank);
  scatter( &gfield, &lfield, 1, root );
}

template<typename DATA_TYPE>
void GatherScatter::pack_send_buffer( const mpl::Field<DATA_TYPE const>& field,
                                      const std::vector<int>& sendmap,
                                      DATA_TYPE send_buffer[] ) const
{
  const int sendcnt = sendmap.size();

  int ibuf = 0;
  const int send_stride = field.var_strides[0]*field.var_extents[0];

  switch( field.var_rank )
  {
  case 1:
    for( int p=0; p<sendcnt; ++p)
    {
      const int pp = send_stride*sendmap[p];
      for( int i=0; i<field.var_extents[0]; ++i )
      {
        DATA_TYPE tmp =  field.data[pp+i*field.var_strides[0]];
        send_buffer[ibuf++] = tmp;
      }
    }
    break;
  case 2:
    for( int p=0; p<sendcnt; ++p)
    {
      const int pp = send_stride*sendmap[p];
      for( int i=0; i<field.var_extents[0]; ++i )
      {
        const int ii = pp + i*field.var_strides[0];
        for( int j=0; j<field.var_extents[1]; ++j )
        {
          send_buffer[ibuf++] = field.data[ii+j*field.var_strides[1]];
        }
      }
    }
    break;
  case 3:
    for( int p=0; p<sendcnt; ++p)
    {
      const int pp = send_stride*sendmap[p];
      for( int i=0; i<field.var_extents[0]; ++i )
      {
        const int ii = pp + i*field.var_strides[0];
        for( int j=0; j<field.var_extents[1]; ++j )
        {
          const int jj = ii + j*field.var_strides[1];
          for( int k=0; k<field.var_extents[2]; ++k )
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
                                        const mpl::Field<DATA_TYPE>& field ) const
{
  const int recvcnt = recvmap.size();
  bool field_changed = false;
  DATA_TYPE tmp;
  int ibuf = 0;
  const int recv_stride = field.var_strides[0]*field.var_extents[0];

  switch( field.var_rank )
  {
  case 1:
    for( int p=0; p<recvcnt; ++p)
    {
      const int pp = recv_stride*recvmap[p];
      for( int i=0; i<field.var_extents[0]; ++i)
      {
        field.data[ pp + i*field.var_strides[0] ] = recv_buffer[ibuf++];
      }
    }
    break;
  case 2:
    for( int p=0; p<recvcnt; ++p)
    {
      const int pp = recv_stride*recvmap[p];
      for( int i=0; i<field.var_extents[0]; ++i )
      {
        const int ii = pp + i*field.var_strides[0];
        for( int j=0; j<field.var_extents[1]; ++j )
        {
          field.data[ ii + j*field.var_strides[1] ] = recv_buffer[ibuf++];
        }
      }
    }
    break;
  case 3:
    for( int p=0; p<recvcnt; ++p)
    {
      const int pp = recv_stride*recvmap[p];
      for( int i=0; i<field.var_extents[0]; ++i )
      {
        const int ii = pp + i*field.var_strides[0];
        for( int j=0; j<field.var_extents[1]; ++j )
        {
          const int jj = ii + j*field.var_strides[1];
          for( int k=0; k<field.var_extents[2]; ++k )
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


//template <typename DATA_TYPE>
//void GatherScatter::gather( const DATA_TYPE ldata[],
//                            DATA_TYPE gdata[],
//                            const int nb_vars ) const
//{
//  int strides[] = {1};
//  int extents[] = {nb_vars};
//  gather( ldata, strides, extents, 1,
//           gdata, strides, extents, 1 );
//}


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
void GatherScatter::gather( const ArrayView<DATA_TYPE,LRANK>& ldata,
                            ArrayView<DATA_TYPE,GRANK>& gdata,
                            const int root ) const
{
  if( ldata.size() == parsize_ && gdata.size() == glbcnt_ )
  {
    std::vector< mpl::Field<DATA_TYPE const> > lfields(1, mpl::Field<DATA_TYPE const>(ldata) );
    std::vector< mpl::Field<DATA_TYPE> >       gfields(1, mpl::Field<DATA_TYPE>(gdata) );
    gather( lfields.data(), gfields.data(), 1, root );
  }
  else
  {
    DEBUG_VAR(parsize_);
    DEBUG_VAR(ldata.size());
    DEBUG_VAR(gdata.size());
    NOTIMP; // Need to implement with parallel ranks > 1
  }
}

template <typename DATA_TYPE, int GRANK, int LRANK>
void GatherScatter::scatter( const ArrayView<DATA_TYPE,GRANK>& gdata,
                             ArrayView<DATA_TYPE,LRANK>& ldata,
                             const int root ) const
{
  if( ldata.size() == parsize_ && gdata.size() == glbcnt_ )
  {
    std::vector< mpl::Field<DATA_TYPE const> > gfields(1, mpl::Field<DATA_TYPE const>(gdata) );
    std::vector< mpl::Field<DATA_TYPE> >       lfields(1, mpl::Field<DATA_TYPE>(ldata) );
    scatter( gfields.data(), lfields.data(), 1, root );
  }
  else
  {
    DEBUG_VAR(ldata.size());
    DEBUG_VAR(gdata.size());
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
  void atlas__GatherScatter__gather_int (GatherScatter* This, int ldata[], int lvar_strides[], int lvar_extents[], int lvar_rank, int gdata[], int gvar_strides[], int gvar_extents[], int gvar_rank);
  void atlas__GatherScatter__gather_float (GatherScatter* This, float ldata[], int lvar_strides[], int lvar_extents[], int lvar_rank, float gdata[], int gvar_strides[], int gvar_extents[], int gvar_rank);
  void atlas__GatherScatter__gather_double (GatherScatter* This, double ldata[], int lvar_strides[], int lvar_extents[], int lvar_rank, double gdata[], int gvar_strides[], int gvar_extents[], int gvar_rank);
  void atlas__GatherScatter__scatter_int (GatherScatter* This, int gdata[], int gvar_strides[], int gvar_extents[], int gvar_rank, int ldata[], int lvar_strides[], int lvar_extents[], int lvar_rank);
  void atlas__GatherScatter__scatter_float (GatherScatter* This, float gdata[], int gvar_strides[], int gvar_extents[], int gvar_rank, float ldata[], int lvar_strides[], int lvar_extents[], int lvar_rank);
  void atlas__GatherScatter__scatter_double (GatherScatter* This, double gdata[], int gvar_strides[], int gvar_extents[], int gvar_rank, double ldata[], int lvar_strides[], int lvar_extents[], int lvar_rank);
}

// ------------------------------------------------------------------

//typedef GatherScatter Gather;

} // namespace mpl
} // namespace atlas

#endif // Gather_h
