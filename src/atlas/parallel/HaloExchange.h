/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

/// @file   HaloExchange.h
/// @author Willem Deconinck
/// @date   Nov 2013

#pragma once

#include <stdexcept>
#include <string>
#include <vector>

#include "atlas/parallel/HaloExchangeImpl.h"
#include "atlas/parallel/mpi/Statistics.h"
#include "atlas/parallel/mpi/mpi.h"


#include "atlas/array/ArrayView.h"
#include "atlas/array/ArrayViewDefs.h"
#include "atlas/array/ArrayViewUtil.h"
#include "atlas/array/SVector.h"
#include "atlas/array_fwd.h"
#include "atlas/library/config.h"
#include "atlas/runtime/Exception.h"
#include "atlas/util/Object.h"

#ifdef ATLAS_GRIDTOOLS_STORAGE_BACKEND_CUDA
#include "atlas/parallel/HaloExchangeCUDA.h"
#endif

namespace atlas {
namespace parallel {

class HaloExchange : public util::Object {
public:
    HaloExchange();
    HaloExchange( const std::string& name );
    virtual ~HaloExchange();

public:  // methods
    const std::string& name() const { return name_; }

    void setup( const int part[], const idx_t remote_idx[], const int base, idx_t size );

    void setup( const int part[], const idx_t remote_idx[], const int base, idx_t size, idx_t halo_begin );

    //  template <typename DATA_TYPE>
    //  void execute( DATA_TYPE field[], idx_t nb_vars ) const;

    template <typename DATA_TYPE, int RANK, typename ParallelDim = array::FirstDim>
    void execute( array::Array& field, bool on_device = false ) const;

private:  // methods
    void create_mappings( std::vector<int>& send_map, std::vector<int>& recv_map, idx_t nb_vars ) const;

    template <int N, int P>
    void create_mappings_impl( std::vector<int>& send_map, std::vector<int>& recv_map, idx_t nb_vars ) const;

    idx_t index( idx_t i, idx_t j, idx_t k, idx_t ni, idx_t nj, idx_t nk ) const { return ( i + ni * ( j + nj * k ) ); }

    idx_t index( idx_t i, idx_t j, idx_t ni, idx_t nj ) const { return ( i + ni * j ); }

    template <int ParallelDim, typename DATA_TYPE, int RANK>
    void pack_send_buffer( const array::ArrayView<DATA_TYPE, RANK>& hfield,
                           const array::ArrayView<DATA_TYPE, RANK>& dfield, DATA_TYPE* send_buffer,
                           int send_buffer_size, const bool on_device ) const;

    template <int ParallelDim, typename DATA_TYPE, int RANK>
    void unpack_recv_buffer( const DATA_TYPE* recv_buffer, int recv_buffer_size,
                             array::ArrayView<DATA_TYPE, RANK>& hfield, array::ArrayView<DATA_TYPE, RANK>& dfield,
                             const bool on_device ) const;

    template <typename DATA_TYPE, int RANK>
    void var_info( const array::ArrayView<DATA_TYPE, RANK>& arr, std::vector<idx_t>& varstrides,
                   std::vector<idx_t>& varshape ) const;

private:  // data
    std::string name_;
    bool is_setup_;

    int sendcnt_;
    int recvcnt_;
    std::vector<int> sendcounts_;
    std::vector<int> senddispls_;
    std::vector<int> recvcounts_;
    std::vector<int> recvdispls_;
    array::SVector<int> sendmap_;
    array::SVector<int> recvmap_;
    int parsize_;

    int nproc;
    int myproc;

public:
    struct Backdoor {
        int parsize;
    } backdoor;
};

template <typename DATA_TYPE, int RANK, typename ParallelDim>
void HaloExchange::execute( array::Array& field, bool on_device ) const {
    if ( !is_setup_ ) {
        throw_Exception( "HaloExchange was not setup", Here() );
    }

    ATLAS_TRACE( "HaloExchange", {"halo-exchange"} );

    auto field_hv = array::make_host_view<DATA_TYPE, RANK>( field );

    int tag                   = 1;
    constexpr int parallelDim = array::get_parallel_dim<ParallelDim>( field_hv );
    idx_t var_size            = array::get_var_size<parallelDim>( field_hv );
    int send_size             = sendcnt_ * var_size;
    int recv_size             = recvcnt_ * var_size;

    DATA_TYPE* send_buffer{nullptr};
    DATA_TYPE* recv_buffer{nullptr};
    if ( on_device ) {
        util::allocate_devicemem( send_buffer, send_size );
        util::allocate_devicemem( recv_buffer, recv_size );
    }
    else {
        util::allocate_hostmem( send_buffer, send_size );
        util::allocate_hostmem( recv_buffer, recv_size );
    }
    std::vector<int> send_displs( nproc );
    std::vector<int> recv_displs( nproc );
    std::vector<int> send_counts( nproc );
    std::vector<int> recv_counts( nproc );

    std::vector<eckit::mpi::Request> send_req( nproc );
    std::vector<eckit::mpi::Request> recv_req( nproc );

    for ( int jproc = 0; jproc < nproc; ++jproc ) {
        send_counts[jproc] = sendcounts_[jproc] * var_size;
        recv_counts[jproc] = recvcounts_[jproc] * var_size;
        send_displs[jproc] = senddispls_[jproc] * var_size;
        recv_displs[jproc] = recvdispls_[jproc] * var_size;
    }

    auto field_dv =
        on_device ? array::make_device_view<DATA_TYPE, RANK>( field ) : array::make_host_view<DATA_TYPE, RANK>( field );

    ATLAS_TRACE_MPI( IRECEIVE ) {
        /// Let MPI know what we like to receive
        for ( int jproc = 0; jproc < nproc; ++jproc ) {
            if ( recv_counts[jproc] > 0 ) {
                recv_req[jproc] =
                    mpi::comm().iReceive( &recv_buffer[recv_displs[jproc]], recv_counts[jproc], jproc, tag );
            }
        }
    }

    /// Pack
    pack_send_buffer<parallelDim>( field_hv, field_dv, send_buffer, send_size, on_device );

    /// Send
    ATLAS_TRACE_MPI( ISEND ) {
        for ( int jproc = 0; jproc < nproc; ++jproc ) {
            if ( send_counts[jproc] > 0 ) {
                send_req[jproc] = mpi::comm().iSend( &send_buffer[send_displs[jproc]], send_counts[jproc], jproc, tag );
            }
        }
    }

    /// Wait for receiving to finish
    ATLAS_TRACE_MPI( WAIT, "mpi-wait receive" ) {
        for ( int jproc = 0; jproc < nproc; ++jproc ) {
            if ( recvcounts_[jproc] > 0 ) {
                mpi::comm().wait( recv_req[jproc] );
            }
        }
    }

    /// Unpack
    unpack_recv_buffer<parallelDim>( recv_buffer, recv_size, field_hv, field_dv, on_device );

    /// Wait for sending to finish
    ATLAS_TRACE_MPI( WAIT, "mpi-wait send" ) {
        for ( int jproc = 0; jproc < nproc; ++jproc ) {
            if ( sendcounts_[jproc] > 0 ) {
                mpi::comm().wait( send_req[jproc] );
            }
        }
    }
    if ( on_device ) {
        util::delete_devicemem( send_buffer );
        util::delete_devicemem( recv_buffer );
    }
    else {
        util::delete_hostmem( send_buffer );
        util::delete_hostmem( recv_buffer );
    }
}

template <int ParallelDim, int RANK>
struct halo_packer {
    template <typename DATA_TYPE>
    static void pack( const int sendcnt, array::SVector<int> const& sendmap,
                      const array::ArrayView<DATA_TYPE, RANK>& field, DATA_TYPE* send_buffer, int send_buffer_size ) {
        idx_t ibuf = 0;
        for ( int node_cnt = 0; node_cnt < sendcnt; ++node_cnt ) {
            const idx_t node_idx = sendmap[node_cnt];
            halo_packer_impl<ParallelDim, RANK, 0>::apply( ibuf, node_idx, field, send_buffer );
        }
    }

    template <typename DATA_TYPE>
    static void unpack( const int recvcnt, array::SVector<int> const& recvmap, const DATA_TYPE* recv_buffer,
                        int recv_buffer_size, array::ArrayView<DATA_TYPE, RANK>& field ) {
        idx_t ibuf = 0;
        for ( int node_cnt = 0; node_cnt < recvcnt; ++node_cnt ) {
            const idx_t node_idx = recvmap[node_cnt];
            halo_unpacker_impl<ParallelDim, RANK, 0>::apply( ibuf, node_idx, recv_buffer, field );
        }
    }
};

template <int ParallelDim, typename DATA_TYPE, int RANK>
void HaloExchange::pack_send_buffer( const array::ArrayView<DATA_TYPE, RANK>& hfield,
                                     const array::ArrayView<DATA_TYPE, RANK>& dfield, DATA_TYPE* send_buffer,
                                     int send_size, const bool on_device ) const {
    ATLAS_TRACE();
#if ATLAS_GRIDTOOLS_STORAGE_BACKEND_CUDA
    if ( on_device ) {
        halo_packer_cuda<ParallelDim, DATA_TYPE, RANK>::pack( sendcnt_, sendmap_, hfield, dfield, send_buffer,
                                                              send_size );
    }
    else
#endif
        halo_packer<ParallelDim, RANK>::pack( sendcnt_, sendmap_, dfield, send_buffer, send_size );
}

template <int ParallelDim, typename DATA_TYPE, int RANK>
void HaloExchange::unpack_recv_buffer( const DATA_TYPE* recv_buffer, int recv_size,
                                       array::ArrayView<DATA_TYPE, RANK>& hfield,
                                       array::ArrayView<DATA_TYPE, RANK>& dfield, const bool on_device ) const {
    ATLAS_TRACE();
#if ATLAS_GRIDTOOLS_STORAGE_BACKEND_CUDA
    if ( on_device ) {
        halo_packer_cuda<ParallelDim, DATA_TYPE, RANK>::unpack( recvcnt_, recvmap_, recv_buffer, recv_size, hfield,
                                                                dfield );
    }
    else
#endif
        halo_packer<ParallelDim, RANK>::unpack( recvcnt_, recvmap_, recv_buffer, recv_size, dfield );
}

// template<typename DATA_TYPE>
// void HaloExchange::execute( DATA_TYPE field[], idx_t nb_vars ) const
//{
//    throw_AssertionFailed("Call not supported");

//  idx_t strides[] = {1};
//  idx_t shape[] = {nb_vars};
//  execute( field, strides, shape, 1);
//}

// template <typename DATA_TYPE, int RANK>
// void HaloExchange::execute( array::ArrayView<DATA_TYPE,RANK>&& field ) const
//{
//    execute(field);
//}

//----------------------------------------------------------------------------------------------------------------------
// C wrapper interfaces to C++ routines
extern "C" {
HaloExchange* atlas__HaloExchange__new();
void atlas__HaloExchange__delete( HaloExchange* This );
void atlas__HaloExchange__setup( HaloExchange* This, int part[], idx_t remote_idx[], int base, int size );
void atlas__HaloExchange__execute_strided_int( HaloExchange* This, int field[], int var_strides[], int var_shape[],
                                               int var_rank );
void atlas__HaloExchange__execute_strided_long( HaloExchange* This, long field[], int var_strides[], int var_shape[],
                                                int var_rank );
void atlas__HaloExchange__execute_strided_float( HaloExchange* This, float field[], int var_strides[], int var_shape[],
                                                 int var_rank );
void atlas__HaloExchange__execute_strided_double( HaloExchange* This, double field[], int var_strides[],
                                                  int var_shape[], int var_rank );
void atlas__HaloExchange__execute_int( HaloExchange* This, int field[], int var_rank );
void atlas__HaloExchange__execute_float( HaloExchange* This, float field[], int var_rank );
void atlas__HaloExchange__execute_double( HaloExchange* This, double field[], int var_rank );
}

//----------------------------------------------------------------------------------------------------------------------

}  // namespace parallel
}  // namespace atlas
