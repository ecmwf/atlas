/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

/// @author Willem Deconinck
/// @date   Nov 2013

#ifndef HaloExchange_h
#define HaloExchange_h

#include <string>
#include <vector>
#include <stdexcept>

#include "eckit/memory/SharedPtr.h"
#include "eckit/memory/Owned.h"
#include "eckit/exception/Exceptions.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/parallel/mpi/Statistics.h"

#include "atlas/array_fwd.h"
#include "atlas/array/ArrayView.h"
#include "atlas/array/SVector.h"
#include "atlas/runtime/Log.h"
#include "atlas/array/ArrayViewDefs.h"
#include "atlas/array/ArrayViewUtil.h"

#ifdef ATLAS_GRIDTOOLS_STORAGE_BACKEND_CUDA
#include "atlas/parallel/HaloExchangeCUDA.h"
#endif


namespace atlas {
namespace parallel {

class HaloExchange: public eckit::Owned {

public: // types
    typedef eckit::SharedPtr<HaloExchange> Ptr;
public:
  HaloExchange();
  HaloExchange(const std::string& name);
  virtual ~HaloExchange() {}

public: // methods

  const std::string& name() const { return name_; }

  void setup( const int part[],
              const int remote_idx[], const int base,
              size_t size );

//  template <typename DATA_TYPE>
//  void execute( DATA_TYPE field[], size_t nb_vars ) const;

  template<typename DATA_TYPE, int RANK, typename ParallelDim = array::FirstDim >
  void execute( array::Array& field, bool on_device = false) const;

private: // methods

  void create_mappings( std::vector<int>& send_map, std::vector<int>& recv_map, size_t nb_vars ) const;

  template<int N, int P>
  void create_mappings_impl( std::vector<int>& send_map, std::vector<int>& recv_map, size_t nb_vars ) const;

  size_t index(size_t i, size_t j, size_t k, size_t ni, size_t nj, size_t nk) const { return( i + ni*( j + nj*k) ); }

  size_t index(size_t i, size_t j, size_t ni, size_t nj) const { return( i + ni*j ); }


  template< int ParallelDim, typename DATA_TYPE, int RANK>
  void pack_send_buffer( const array::ArrayView<DATA_TYPE, RANK, array::Intent::ReadOnly>& hfield,
                         const array::ArrayView<DATA_TYPE, RANK>& dfield,
                         array::SVector<DATA_TYPE>& send_buffer ) const;

  template< int ParallelDim, typename DATA_TYPE, int RANK>
  void unpack_recv_buffer(const array::SVector<DATA_TYPE>& recv_buffer,
                          const array::ArrayView<DATA_TYPE, RANK,array::Intent::ReadOnly>& hfield,
                          array::ArrayView<DATA_TYPE, RANK>& dfield) const;

  template<typename DATA_TYPE, int RANK>
  void var_info( const array::ArrayView<DATA_TYPE,RANK>& arr,
                 std::vector<size_t>& varstrides,
                 std::vector<size_t>& varshape ) const;

private: // data
  std::string name_;
  bool is_setup_;

  int               sendcnt_;
  int               recvcnt_;
  std::vector<int>  sendcounts_;
  std::vector<int>  senddispls_;
  std::vector<int>  recvcounts_;
  std::vector<int>  recvdispls_;
  array::SVector<int>  sendmap_;
  array::SVector<int>  recvmap_;
  int parsize_;

  int nproc;
  int myproc;

  std::vector<int> bounds_;
  int par_bound_;

};

template<typename DATA_TYPE, int RANK, typename ParallelDim>
void HaloExchange::execute(array::Array& field, bool on_device) const
{
  if( ! is_setup_ )
  {
    throw eckit::SeriousBug("HaloExchange was not setup",Here());
  }

  ATLAS_TRACE("HaloExchange",{"halo-exchange"});

  auto field_hv = array::make_host_view<DATA_TYPE, RANK, array::Intent::ReadOnly>(field);

  int tag=1;
  constexpr int parallelDim = array::get_parallel_dim<ParallelDim>(field_hv);
  size_t var_size = array::get_var_size< parallelDim >(field_hv);
  int send_size  = sendcnt_ * var_size;
  int recv_size  = recvcnt_ * var_size;

  array::SVector<DATA_TYPE  > send_buffer(send_size);
  array::SVector<DATA_TYPE  > recv_buffer(recv_size);
  std::vector<int        > send_displs(nproc    );
  std::vector<int        > recv_displs(nproc    );
  std::vector<int        > send_counts(nproc    );
  std::vector<int        > recv_counts(nproc    );

  std::vector<eckit::mpi::Request> send_req(nproc    );
  std::vector<eckit::mpi::Request> recv_req(nproc    );

  for (int jproc=0; jproc < nproc; ++jproc)
  {
    send_counts[jproc] = sendcounts_[jproc]*var_size;
    recv_counts[jproc] = recvcounts_[jproc]*var_size;
    send_displs[jproc] = senddispls_[jproc]*var_size;
    recv_displs[jproc] = recvdispls_[jproc]*var_size;
  }

  auto field_dv = on_device ? array::make_device_view<DATA_TYPE, RANK>(field) :
      array::make_host_view<DATA_TYPE, RANK>(field);

  ATLAS_TRACE_MPI( IRECEIVE ) {
    /// Let MPI know what we like to receive
    for(int jproc=0; jproc < nproc; ++jproc)
    {
      if(recv_counts[jproc] > 0)
      {
          recv_req[jproc] = parallel::mpi::comm().iReceive(&recv_buffer[recv_displs[jproc]], recv_counts[jproc], jproc, tag);
      }
    }
  }

  /// Pack
  pack_send_buffer<parallelDim>(field_hv, field_dv,send_buffer);

  /// Send
  ATLAS_TRACE_MPI( ISEND ) {
    for(int jproc=0; jproc < nproc; ++jproc)
    {
      if(send_counts[jproc] > 0)
      {
         send_req[jproc] = parallel::mpi::comm().iSend(
              &send_buffer[send_displs[jproc]],
              send_counts[jproc], jproc, tag);
      }
    }
  }

  /// Wait for receiving to finish
  ATLAS_TRACE_MPI( WAIT, "mpi-wait receive" ) {
    for (int jproc=0; jproc < nproc; ++jproc)
    {
      if( recvcounts_[jproc] > 0)
      {
          parallel::mpi::comm().wait(recv_req[jproc]);
      }
    }
  }

  /// Unpack
  unpack_recv_buffer<parallelDim>(recv_buffer, field_hv, field_dv);

  /// Wait for sending to finish
  ATLAS_TRACE_MPI( WAIT, "mpi-wait send" ) {
    for (int jproc=0; jproc < nproc; ++jproc)
    {
      if( sendcounts_[jproc] > 0)
      {
          parallel::mpi::comm().wait(send_req[jproc]);
      }
    }
  }

}

template<int ParallelDim, int Cnt, int CurrentDim>
struct halo_packer_impl {
    template<typename DATA_TYPE, int RANK, typename ... Idx>
    static void apply(size_t& buf_idx, const size_t node_idx, const array::ArrayView<DATA_TYPE, RANK, array::Intent::ReadOnly>& field,
                      array::SVector<DATA_TYPE>& send_buffer, Idx... idxs) {
      for( size_t i=0; i< field.template shape<CurrentDim>(); ++i ) {
        halo_packer_impl<ParallelDim, Cnt-1, CurrentDim+1>::apply(buf_idx, node_idx, field, send_buffer, idxs..., i);
      }
    }
};

template<int ParallelDim, int Cnt>
struct halo_packer_impl<ParallelDim, Cnt, ParallelDim>
{
    template<typename DATA_TYPE, int RANK, typename ... Idx>
    static void apply(size_t& buf_idx, const size_t node_idx, const array::ArrayView<DATA_TYPE, RANK, array::Intent::ReadOnly>& field,
                      array::SVector<DATA_TYPE>& send_buffer, Idx... idxs) {
        //TODO node should be inserted in the right place accordim to parallel dim
      halo_packer_impl<ParallelDim, Cnt-1, ParallelDim+1>::apply(buf_idx, node_idx, field, send_buffer, idxs...);
    }
};

template<int ParallelDim, int CurrentDim>
struct halo_packer_impl<ParallelDim, 0, CurrentDim> {
    template<typename DATA_TYPE, int RANK, typename ...Idx>
    static void apply(size_t& buf_idx, size_t node_idx, const array::ArrayView<DATA_TYPE, RANK, array::Intent::ReadOnly>& field,
                     array::SVector<DATA_TYPE>& send_buffer, Idx...idxs)
    {
      send_buffer[buf_idx++] = field(node_idx, idxs...);
    }
};

template<int ParallelDim, int Cnt, int CurrentDim>
struct halo_unpacker_impl {
    template<typename DATA_TYPE, int RANK, typename ... Idx>
    static void apply(size_t& buf_idx, const size_t node_idx, array::SVector<DATA_TYPE> const & recv_buffer, 
                      array::ArrayView<DATA_TYPE, RANK>& field, Idx... idxs) {
        for( size_t i=0; i< field.template shape<CurrentDim>(); ++i ) {
            halo_unpacker_impl<ParallelDim, Cnt-1, CurrentDim+1>::apply(buf_idx, node_idx, recv_buffer, field, idxs..., i);
        }
    }
};

template<int ParallelDim, int Cnt>
struct halo_unpacker_impl<ParallelDim, Cnt, ParallelDim> {
    template<typename DATA_TYPE, int RANK, typename ... Idx>
    static void apply(size_t& buf_idx, const size_t node_idx, array::SVector<DATA_TYPE> const & recv_buffer,
                      array::ArrayView<DATA_TYPE, RANK>& field, Idx... idxs) {
        halo_unpacker_impl<ParallelDim, Cnt-1, ParallelDim+1>::apply(buf_idx, node_idx, recv_buffer, field, idxs...);
    }
};

template<int ParallelDim, int CurrentDim>
struct halo_unpacker_impl<ParallelDim, 0, CurrentDim> {
    template<typename DATA_TYPE, int RANK, typename ...Idx>
    static void apply(size_t& buf_idx, size_t node_idx, array::SVector<DATA_TYPE> const & recv_buffer,
                     array::ArrayView<DATA_TYPE, RANK>& field, Idx...idxs)
    {
      field(node_idx, idxs...) = recv_buffer[buf_idx++];
    }
};

template<int ParallelDim, int RANK>
struct halo_packer {
    template<typename DATA_TYPE>
    static void pack(const unsigned int sendcnt, array::SVector<int> const & sendmap,
                     const array::ArrayView<DATA_TYPE, RANK, array::Intent::ReadOnly>& field, array::SVector<DATA_TYPE>& send_buffer )
    {
      size_t ibuf = 0;
      for(int node_cnt=0; node_cnt < sendcnt; ++node_cnt)
      {
        const size_t node_idx = sendmap[node_cnt];
        halo_packer_impl<ParallelDim, RANK,0>::apply(ibuf, node_idx, field, send_buffer);
      }
    }

    template<typename DATA_TYPE>
    static void unpack(const unsigned int recvcnt, array::SVector<int> const & recvmap,
                     array::SVector<DATA_TYPE> const & recv_buffer, array::ArrayView<DATA_TYPE, RANK>& field )
    {
      size_t ibuf = 0;
      for(int node_cnt=0; node_cnt < recvcnt; ++node_cnt)
      {
        const size_t node_idx = recvmap[node_cnt];
        halo_unpacker_impl<ParallelDim, RANK,0>::apply(ibuf, node_idx, recv_buffer, field);
      }
    }

};

template<int ParallelDim, typename DATA_TYPE, int RANK>
void HaloExchange::pack_send_buffer( const array::ArrayView<DATA_TYPE, RANK, array::Intent::ReadOnly>& hfield,
                                     const array::ArrayView<DATA_TYPE, RANK>& dfield,
                                     array::SVector<DATA_TYPE>& send_buffer ) const
{
  ATLAS_TRACE();
#if ATLAS_GRIDTOOLS_STORAGE_BACKEND_CUDA
    halo_packer_cuda<DATA_TYPE, RANK>::pack(sendcnt_, sendmap_, hfield, dfield, send_buffer);
#else
    halo_packer<ParallelDim, RANK>::pack(sendcnt_, sendmap_, hfield, send_buffer);
#endif
}

template<int ParallelDim, typename DATA_TYPE, int RANK>
void HaloExchange::unpack_recv_buffer( const array::SVector<DATA_TYPE>& recv_buffer,
                                       const array::ArrayView<DATA_TYPE, RANK, array::Intent::ReadOnly>& hfield,
                                       array::ArrayView<DATA_TYPE,RANK>& dfield) const
{
  ATLAS_TRACE();

#if ATLAS_GRIDTOOLS_STORAGE_BACKEND_CUDA
    halo_packer_cuda<DATA_TYPE, RANK>::unpack(recvcnt_, recvmap_, recv_buffer, hfield, dfield);
#else
    halo_packer<ParallelDim, RANK>::unpack(recvcnt_, recvmap_, recv_buffer, dfield);
#endif

}

//template<typename DATA_TYPE>
//void HaloExchange::execute( DATA_TYPE field[], size_t nb_vars ) const
//{
//    throw eckit::AssertionFailed("Call not supported");

//  size_t strides[] = {1};
//  size_t shape[] = {nb_vars};
//  execute( field, strides, shape, 1);
//}


//template <typename DATA_TYPE, int RANK>
//void HaloExchange::execute( array::ArrayView<DATA_TYPE,RANK>&& field ) const
//{
//    execute(field);
//}

//----------------------------------------------------------------------------------------------------------------------
// C wrapper interfaces to C++ routines
extern "C"
{
  HaloExchange* atlas__HaloExchange__new ();
  void atlas__HaloExchange__delete (HaloExchange* This);
  void atlas__HaloExchange__setup (HaloExchange* This, int part[], int remote_idx[], int base, int size);
  void atlas__HaloExchange__execute_strided_int (HaloExchange* This, int field[], int var_strides[], int var_shape[], int var_rank);
  void atlas__HaloExchange__execute_strided_long (HaloExchange* This, long field[], int var_strides[], int var_shape[], int var_rank);
  void atlas__HaloExchange__execute_strided_float (HaloExchange* This, float field[], int var_strides[], int var_shape[], int var_rank);
  void atlas__HaloExchange__execute_strided_double (HaloExchange* This, double field[], int var_strides[], int var_shape[], int var_rank);
  void atlas__HaloExchange__execute_int (HaloExchange* This, int field[], int var_rank);
  void atlas__HaloExchange__execute_float (HaloExchange* This, float field[], int var_rank);
  void atlas__HaloExchange__execute_double (HaloExchange* This, double field[], int var_rank);

}

//----------------------------------------------------------------------------------------------------------------------

} // namespace parallel
} // namespace atlas

#endif // HaloExchange_h
