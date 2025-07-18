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

#include "atlas/parallel/detail/Packer.h"
#include "atlas/parallel/mpi/Statistics.h"
#include "atlas/parallel/mpi/mpi.h"

#include "atlas/array/ArrayView.h"
#include "atlas/array/ArrayViewDefs.h"
#include "atlas/array/ArrayViewUtil.h"
#include "atlas/array/SVector.h"
#include "atlas/array_fwd.h"
#include "atlas/library/config.h"
#include "atlas/library/defines.h"
#include "atlas/runtime/Exception.h"
#include "atlas/util/Object.h"

namespace atlas {
namespace parallel {

class HaloExchange : public util::Object {
public:
    HaloExchange();
    HaloExchange(const std::string& name);

    virtual ~HaloExchange();

public:  // methods
    const std::string& name() const { return name_; }

    void setup(const int part[], const idx_t remote_idx[], const int base, idx_t size);
    void setup(const std::string& mpi_comm, const int part[], const idx_t remote_idx[], const int base, idx_t size);

    void setup(const int part[], const idx_t remote_idx[], const int base, idx_t size, idx_t halo_begin);
    void setup(const std::string& mpi_comm, const int part[], const idx_t remote_idx[], const int base, idx_t size, idx_t halo_begin);

    template <typename DATA_TYPE, int RANK, typename ParallelDim = array::FirstDim>
    void execute(array::Array& field, bool on_device = false) const;

    template <typename DATA_TYPE, int RANK, typename ParallelDim = array::FirstDim>
    void execute_adjoint(array::Array& field, bool on_device = false) const;

private:  // methods
    idx_t index(idx_t i, idx_t j, idx_t k, idx_t ni, idx_t nj, idx_t /*nk*/) const { return (i + ni * (j + nj * k)); }

    idx_t index(idx_t i, idx_t j, idx_t ni, idx_t /*nj*/) const { return (i + ni * j); }

    template <typename DATA_TYPE>
    void counts_displs_setup(const idx_t var_size, std::vector<int>& send_counts_init,
                             std::vector<int>& recv_counts_init, std::vector<int>& send_counts,
                             std::vector<int>& recv_counts, std::vector<int>& send_displs,
                             std::vector<int>& recv_displs) const;


    template <typename DATA_TYPE>
    void ireceive(int tag, std::vector<int>& recv_displs, std::vector<int>& recv_counts,
                  std::vector<eckit::mpi::Request>& recv_req, DATA_TYPE* recv_buffer) const;

    template <typename DATA_TYPE>
    void isend_and_wait_for_receive(int tag, std::vector<int>& recv_counts_init,
                                    std::vector<eckit::mpi::Request>& recv_req, std::vector<int>& send_displs,
                                    std::vector<int>& send_counts, std::vector<eckit::mpi::Request>& send_req,
                                    DATA_TYPE* send_buffer) const;

    void wait_for_send(std::vector<int>& send_counts, std::vector<eckit::mpi::Request>& send_req) const;

    template <typename DATA_TYPE>
    DATA_TYPE* allocate_buffer(const int buffer_size, const bool on_device) const;

    template <typename DATA_TYPE>
    void deallocate_buffer(DATA_TYPE* buffer, const int buffer_size, const bool on_device) const;

    template <int ParallelDim, typename DATA_TYPE, int RANK>
    void pack_send_buffer(const array::ArrayView<DATA_TYPE, RANK>& array,
                          DATA_TYPE* send_buffer, int send_buffer_size,
                          const bool on_device) const;

    template <int ParallelDim, typename DATA_TYPE, int RANK>
    void unpack_recv_buffer(const DATA_TYPE* recv_buffer, int recv_buffer_size,
                            array::ArrayView<DATA_TYPE, RANK>& array,
                            const bool on_device) const;

    template <int ParallelDim, typename DATA_TYPE, int RANK>
    void pack_recv_adjoint_buffer(const array::ArrayView<DATA_TYPE, RANK>& array,
                                  DATA_TYPE* recv_buffer, int recv_buffer_size,
                                  const bool on_device) const;

    template <int ParallelDim, typename DATA_TYPE, int RANK>
    void unpack_send_adjoint_buffer(const DATA_TYPE* send_buffer, int send_buffer_size,
                                    array::ArrayView<DATA_TYPE, RANK>& array,
                                    const bool on_device) const;

    template <int ParallelDim, typename DATA_TYPE, int RANK>
    void zero_halos(array::ArrayView<DATA_TYPE, RANK>& array,
                    DATA_TYPE* recv_buffer, int recv_buffer_size,
                    const bool on_device) const;

    template <typename DATA_TYPE, int RANK>
    void var_info(const array::ArrayView<DATA_TYPE, RANK>& arr, std::vector<idx_t>& varstrides,
                  std::vector<idx_t>& varshape) const;

    const mpi::Comm& comm() const {
        return *comm_;
    }

    int devices() const;
    bool is_device_accessible(const void* ptr) const;

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
    const mpi::Comm* comm_;

public:
    struct Backdoor {
        int parsize;
    } backdoor;
};

template <typename DATA_TYPE, int RANK, typename ParallelDim>
void HaloExchange::execute(array::Array& field, bool on_device) const {
    ATLAS_TRACE("HaloExchange", {"halo-exchange"});
    if (!is_setup_) {
        throw_Exception("HaloExchange was not setup", Here());
    }

    on_device = on_device && devices() > 0;

    auto array_view =
        on_device ? array::make_device_view<DATA_TYPE, RANK>(field) : array::make_host_view<DATA_TYPE, RANK>(field);

    constexpr int parallelDim = array::get_parallel_dim<ParallelDim>(array_view);
    idx_t var_size            = array::get_var_size<parallelDim>(array_view);

    int tag(1);
    std::size_t nproc_loc(static_cast<std::size_t>(nproc));
    std::vector<int> inner_counts(nproc_loc), halo_counts(nproc_loc);
    std::vector<int> inner_counts_init(nproc_loc), halo_counts_init(nproc_loc);
    std::vector<int> inner_displs(nproc_loc), halo_displs(nproc_loc);
    std::vector<eckit::mpi::Request> inner_req(nproc_loc), halo_req(nproc_loc);

    int inner_size          = sendcnt_ * var_size;
    int halo_size           = recvcnt_ * var_size;
    DATA_TYPE* inner_buffer = allocate_buffer<DATA_TYPE>(inner_size, on_device);
    DATA_TYPE* halo_buffer  = allocate_buffer<DATA_TYPE>(halo_size, on_device);

    if (on_device) {
        ATLAS_ASSERT( is_device_accessible(inner_buffer) );
        ATLAS_ASSERT( is_device_accessible(halo_buffer) );
        ATLAS_ASSERT( is_device_accessible(array_view.data()) );
    }

    counts_displs_setup<DATA_TYPE>(var_size, inner_counts_init, halo_counts_init, inner_counts, halo_counts,
                                   inner_displs, halo_displs);

    ireceive<DATA_TYPE>(tag, halo_displs, halo_counts, halo_req, halo_buffer);

    /// Pack
    pack_send_buffer<parallelDim>(array_view, inner_buffer, inner_size, on_device);

    isend_and_wait_for_receive<DATA_TYPE>(tag, halo_counts_init, halo_req, inner_displs, inner_counts, inner_req,
                                          inner_buffer);

    /// Unpack
    unpack_recv_buffer<parallelDim>(halo_buffer, halo_size, array_view, on_device);

    wait_for_send(inner_counts_init, inner_req);

    deallocate_buffer<DATA_TYPE>(inner_buffer, inner_size, on_device);
    deallocate_buffer<DATA_TYPE>(halo_buffer, halo_size, on_device);
}

template <typename DATA_TYPE, int RANK, typename ParallelDim>
void HaloExchange::execute_adjoint(array::Array& field, bool on_device) const {
    if (!is_setup_) {
        throw_Exception("HaloExchange was not setup", Here());
    }

    ATLAS_TRACE("HaloExchange", {"halo-exchange-adjoint"});

    on_device = on_device && devices() > 0;

    auto array_view =
        on_device ? array::make_device_view<DATA_TYPE, RANK>(field) : array::make_host_view<DATA_TYPE, RANK>(field);

    constexpr int parallelDim = array::get_parallel_dim<ParallelDim>(array_view);
    idx_t var_size            = array::get_var_size<parallelDim>(array_view);

    int tag(1);
    std::size_t nproc_loc(static_cast<std::size_t>(nproc));
    std::vector<int> halo_counts(nproc_loc), inner_counts(nproc_loc);
    std::vector<int> halo_counts_init(nproc_loc), inner_counts_init(nproc_loc);
    std::vector<int> halo_displs(nproc_loc), inner_displs(nproc_loc);
    std::vector<eckit::mpi::Request> halo_req(nproc_loc), inner_req(nproc_loc);

    int halo_size           = sendcnt_ * var_size;
    int inner_size          = recvcnt_ * var_size;
    DATA_TYPE* halo_buffer  = allocate_buffer<DATA_TYPE>(halo_size, on_device);
    DATA_TYPE* inner_buffer = allocate_buffer<DATA_TYPE>(inner_size, on_device);

    counts_displs_setup<DATA_TYPE>(var_size, halo_counts_init, inner_counts_init, halo_counts, inner_counts,
                                   halo_displs, inner_displs);

    ireceive<DATA_TYPE>(tag, halo_displs, halo_counts, halo_req, halo_buffer);

    /// Pack
    pack_recv_adjoint_buffer<parallelDim>(array_view, inner_buffer, inner_size, on_device);

    /// Send
    isend_and_wait_for_receive<DATA_TYPE>(tag, halo_counts_init, halo_req, inner_displs, inner_counts, inner_req,
                                          inner_buffer);

    /// Unpack
    unpack_send_adjoint_buffer<parallelDim>(halo_buffer, halo_size, array_view, on_device);

    /// Wait for sending to finish
    wait_for_send(inner_counts_init, inner_req);

    zero_halos<parallelDim>(array_view, halo_buffer, halo_size, on_device);

    deallocate_buffer<DATA_TYPE>(halo_buffer, halo_size, on_device);
    deallocate_buffer<DATA_TYPE>(inner_buffer, inner_size, on_device);
}

template <typename DATA_TYPE>
DATA_TYPE* HaloExchange::allocate_buffer(const int buffer_size, const bool on_device) const {
    DATA_TYPE* buffer{nullptr};

    if (on_device) {
        util::allocate_devicemem(buffer, buffer_size);
    }
    else {
        util::allocate_hostmem(buffer, buffer_size);
    }

    return buffer;
}


template <typename DATA_TYPE>
void HaloExchange::deallocate_buffer(DATA_TYPE* buffer, const int buffer_size, const bool on_device) const {
    if (on_device) {
        util::delete_devicemem(buffer, buffer_size);
    }
    else {
        util::delete_hostmem(buffer, buffer_size);
    }
}


template <typename DATA_TYPE>
void HaloExchange::counts_displs_setup(const idx_t var_size, std::vector<int>& send_counts_init,
                                       std::vector<int>& recv_counts_init, std::vector<int>& send_counts,
                                       std::vector<int>& recv_counts, std::vector<int>& send_displs,
                                       std::vector<int>& recv_displs) const {
    for (size_t jproc = 0; jproc < static_cast<size_t>(nproc); ++jproc) {
        send_counts_init[jproc] = sendcounts_[jproc];
        recv_counts_init[jproc] = recvcounts_[jproc];
        send_counts[jproc]      = sendcounts_[jproc] * var_size;
        recv_counts[jproc]      = recvcounts_[jproc] * var_size;
        send_displs[jproc]      = senddispls_[jproc] * var_size;
        recv_displs[jproc]      = recvdispls_[jproc] * var_size;
    }
}

template <typename DATA_TYPE>
void HaloExchange::ireceive(int tag, std::vector<int>& recv_displs, std::vector<int>& recv_counts,
                            std::vector<eckit::mpi::Request>& recv_req, DATA_TYPE* recv_buffer) const {
    ATLAS_TRACE_MPI(IRECEIVE) {
        /// Let MPI know what we like to receive
        for (size_t jproc = 0; jproc < static_cast<size_t>(nproc); ++jproc) {
            if (recv_counts[jproc] > 0) {
                recv_req[jproc] =
                    comm().iReceive(recv_buffer+recv_displs[jproc], recv_counts[jproc], jproc, tag);
            }
        }
    }
}

template <typename DATA_TYPE>
void HaloExchange::isend_and_wait_for_receive(int tag, std::vector<int>& recv_counts_init,
                                              std::vector<eckit::mpi::Request>& recv_req, std::vector<int>& send_displs,
                                              std::vector<int>& send_counts, std::vector<eckit::mpi::Request>& send_req,
                                              DATA_TYPE* send_buffer) const {
    /// Send
    ATLAS_TRACE_MPI(ISEND) {
        for (size_t jproc = 0; jproc < static_cast<size_t>(nproc); ++jproc) {
            if (send_counts[jproc] > 0) {
                send_req[jproc] = comm().iSend(send_buffer+send_displs[jproc], send_counts[jproc], jproc, tag);
            }
        }
    }

    /// Wait for receiving to finish
    ATLAS_TRACE_MPI(WAIT, "mpi-wait receive") {
        for (size_t jproc = 0; jproc < static_cast<size_t>(nproc); ++jproc) {
            if (recv_counts_init[jproc] > 0) {
                comm().wait(recv_req[jproc]);
            }
        }
    }
}

template <int ParallelDim, typename DATA_TYPE, int RANK>
void HaloExchange::zero_halos(array::ArrayView<DATA_TYPE, RANK>& array,
                              DATA_TYPE* recv_buffer, int recv_size,
                              ATLAS_MAYBE_UNUSED const bool on_device) const {
    ATLAS_TRACE();
    detail::Zeroer<ParallelDim, DATA_TYPE, RANK>::zero(recvcnt_, recvmap_, array, recv_buffer, recv_size, on_device);
}

template <int ParallelDim, typename DATA_TYPE, int RANK>
void HaloExchange::pack_send_buffer(const array::ArrayView<DATA_TYPE, RANK>& array,
                                    DATA_TYPE* send_buffer, int send_size,
                                    const bool on_device) const {
    ATLAS_TRACE();
    detail::Packer<ParallelDim, DATA_TYPE, RANK>::pack(sendcnt_, sendmap_, array, send_buffer,
                                                       send_size, on_device);
}

template <int ParallelDim, typename DATA_TYPE, int RANK>
void HaloExchange::unpack_recv_buffer(const DATA_TYPE* recv_buffer, int recv_size,
                                      array::ArrayView<DATA_TYPE, RANK>& array,
                                      const bool on_device) const {
    ATLAS_TRACE();
    detail::Packer<ParallelDim, DATA_TYPE, RANK>::unpack(recvcnt_, recvmap_, recv_buffer, recv_size, array, on_device);
}

template <int ParallelDim, typename DATA_TYPE, int RANK>
void HaloExchange::pack_recv_adjoint_buffer(const array::ArrayView<DATA_TYPE, RANK>& array,
                                            DATA_TYPE* recv_buffer, int recv_size, const bool on_device) const {
    ATLAS_TRACE();
    detail::AdjointPacker<ParallelDim, DATA_TYPE, RANK>::pack(recvcnt_, recvmap_, array, recv_buffer, recv_size, on_device);
}

template <int ParallelDim, typename DATA_TYPE, int RANK>
void HaloExchange::unpack_send_adjoint_buffer(const DATA_TYPE* send_buffer, int send_size,
                                              array::ArrayView<DATA_TYPE, RANK>& array,
                                              ATLAS_MAYBE_UNUSED const bool on_device) const {
    ATLAS_TRACE();
    detail::AdjointPacker<ParallelDim, DATA_TYPE, RANK>::unpack(sendcnt_, sendmap_, send_buffer, send_size, array, on_device);
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
void atlas__HaloExchange__delete(HaloExchange* This);
void atlas__HaloExchange__setup(HaloExchange* This, int part[], idx_t remote_idx[], int base, int size);
void atlas__HaloExchange__execute_strided_int(HaloExchange* This, int field[], int var_strides[], int var_shape[],
                                              int var_rank);
void atlas__HaloExchange__execute_strided_long(HaloExchange* This, long field[], int var_strides[], int var_shape[],
                                               int var_rank);
void atlas__HaloExchange__execute_strided_float(HaloExchange* This, float field[], int var_strides[], int var_shape[],
                                                int var_rank);
void atlas__HaloExchange__execute_strided_double(HaloExchange* This, double field[], int var_strides[], int var_shape[],
                                                 int var_rank);
void atlas__HaloExchange__execute_int(HaloExchange* This, int field[], int var_rank);
void atlas__HaloExchange__execute_float(HaloExchange* This, float field[], int var_rank);
void atlas__HaloExchange__execute_double(HaloExchange* This, double field[], int var_rank);

void atlas__HaloExchange__execute_adjoint_strided_int(HaloExchange* This, int field[], int var_strides[],
                                                      int var_shape[], int var_rank);
void atlas__HaloExchange__execute_adjoint_strided_long(HaloExchange* This, long field[], int var_strides[],
                                                       int var_shape[], int var_rank);
void atlas__HaloExchange__execute_adjoint_strided_float(HaloExchange* This, float field[], int var_strides[],
                                                        int var_shape[], int var_rank);
void atlas__HaloExchange__execute_adjoint_strided_double(HaloExchange* This, double field[], int var_strides[],
                                                         int var_shape[], int var_rank);
void atlas__HaloExchange__execute_adjoint_int(HaloExchange* This, int field[], int var_rank);
void atlas__HaloExchange__execute_adjoint_float(HaloExchange* This, float field[], int var_rank);
void atlas__HaloExchange__execute_adjoint_double(HaloExchange* This, double field[], int var_rank);
}

//----------------------------------------------------------------------------------------------------------------------

}  // namespace parallel
}  // namespace atlas
