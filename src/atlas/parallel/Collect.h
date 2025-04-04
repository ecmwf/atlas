/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>
#include <vector>

#include "atlas/array/Array.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/parallel/HaloExchange.h"

namespace atlas::parallel {

class Collect {
private:
bool is_setup_;

int sendcnt_;
int recvcnt_;
std::vector<int> sendcounts_;
std::vector<int> senddispls_;
std::vector<int> recvcounts_;
std::vector<int> recvdispls_;
array::SVector<int> sendmap_;
array::SVector<int> recvmap_;

int mpi_size_;
int mpi_rank_;
const mpi::Comm* comm_;

public:

const mpi::Comm& comm() const {
    return *comm_;
}

void setup(const idx_t recv_size, const int recv_part[], const idx_t recv_remote_idx[], const int recv_idx_base) {
    setup(mpi::comm().name(), recv_size, recv_part, recv_remote_idx, recv_idx_base );
}

void setup(const std::string& mpi_comm,
    const idx_t recv_size, const int recv_part[], const idx_t recv_remote_idx[], const int recv_idx_base) {
    
    ATLAS_TRACE("Collect::setup()");
    comm_ = &mpi::comm(mpi_comm);
    mpi_size_ = comm().size();
    mpi_rank_ = comm().rank();
    int nproc = mpi_size_;

    recvcnt_ = recv_size;
    sendcounts_.resize(nproc);
    sendcounts_.assign(nproc, 0);
    recvcounts_.resize(nproc);
    recvcounts_.assign(nproc, 0);
    senddispls_.resize(nproc);
    senddispls_.assign(nproc, 0);
    recvdispls_.resize(nproc);
    recvdispls_.assign(nproc, 0);

    atlas_omp_parallel_for(std::size_t jrecv = 0; jrecv < recvcnt_; ++jrecv) {
        int p = recv_part[jrecv];
        atlas_omp_critical {
            ++recvcounts_[p];
        }
    }

    /*
    Find the amount of nodes this rank has to send to each other ranks
    */
    ATLAS_TRACE_MPI(ALLTOALL) { comm().allToAll(recvcounts_, sendcounts_); }

    sendcnt_ = std::accumulate(sendcounts_.begin(), sendcounts_.end(), 0);

    recvdispls_[0] = 0;
    senddispls_[0] = 0;
    for (int jproc = 1; jproc < nproc; ++jproc)  // start at 1
    {
        recvdispls_[jproc] = recvcounts_[jproc - 1] + recvdispls_[jproc - 1];
        senddispls_[jproc] = sendcounts_[jproc - 1] + senddispls_[jproc - 1];
    }

    /*
    Fill vector "send_requests" with remote index of nodes needed, but are on
    other procs
    We can also fill in the vector "recvmap_" which holds local indices of
    requested nodes
    */
    std::vector<int> send_requests(recvcnt_);
    std::vector<int> recv_requests(sendcnt_);
    std::vector<int> cnt(nproc, 0);
    recvmap_.resize(recvcnt_);
#ifdef __PGI
    // No idea why PGI compiler (20.7) in Release build ( -fast -O3 ) decides to vectorize following loop
#pragma loop novector
#endif
    for (int jrecv = 0; jrecv < recvcnt_; ++jrecv) {
        const int p            = recv_part[jrecv];
        const int req_idx      = recvdispls_[p] + cnt[p];
        send_requests[req_idx] = recv_remote_idx[jrecv] - recv_idx_base;
        recvmap_[req_idx]      = jrecv;
        cnt[p]++;
    }

    /*
    Fill vector "recv_requests" with what is needed by other procs
    */
    ATLAS_TRACE_MPI(ALLTOALL) {
        comm().allToAllv(send_requests.data(), recvcounts_.data(), recvdispls_.data(), recv_requests.data(),
                              sendcounts_.data(), senddispls_.data());
    }


    /*
    What needs to be sent to other procs is asked by remote_idx, which is local
    here
    */
    sendmap_.resize(sendcnt_);
    for (int jj = 0; jj < sendcnt_; ++jj) {
        sendmap_[jj] = recv_requests[jj];
    }

    is_setup_        = true;
}

template <typename DATA_TYPE, int RANK, typename ParallelDim = array::FirstDim>
void execute(array::Array& send, array::Array& recv, bool on_device = false) const {
    ATLAS_TRACE("Collect", {"collect"});
    if (!is_setup_) {
        throw_Exception("Collect was not setup", Here());
    }

    if( send.rank() != recv.rank() ){
        throw_Exception("send and recv arrays are of different rank");
    }

    auto send_hv = array::make_host_view<DATA_TYPE, RANK>(send);
    auto recv_hv = array::make_host_view<DATA_TYPE, RANK>(recv);
    auto send_dv =
        on_device ? array::make_device_view<DATA_TYPE, RANK>(send) : array::make_host_view<DATA_TYPE, RANK>(send);
    auto recv_dv =
        on_device ? array::make_device_view<DATA_TYPE, RANK>(recv) : array::make_host_view<DATA_TYPE, RANK>(recv);

    constexpr int parallelDim = array::get_parallel_dim<ParallelDim>(send_hv);
    idx_t var_size            = array::get_var_size<parallelDim>(send_hv);

    int tag(1);
    std::size_t nproc_loc(static_cast<std::size_t>(mpi_size_));
    std::vector<int> send_counts(nproc_loc), recv_counts(nproc_loc);
    std::vector<int> send_counts_init(nproc_loc), recv_counts_init(nproc_loc);
    std::vector<int> send_displs(nproc_loc), recv_displs(nproc_loc);
    std::vector<eckit::mpi::Request> send_req(nproc_loc), recv_req(nproc_loc);

    int send_size          = sendcnt_ * var_size;
    int recv_size          = recvcnt_ * var_size;
    DATA_TYPE* send_buffer = allocate_buffer<DATA_TYPE>(send_size, on_device);
    DATA_TYPE* recv_buffer = allocate_buffer<DATA_TYPE>(recv_size, on_device);

    counts_displs_setup(var_size, send_counts_init, recv_counts_init, send_counts, recv_counts,
                        send_displs, recv_displs);

    ireceive<DATA_TYPE>(tag, recv_displs, recv_counts, recv_req, recv_buffer);

    /// Pack
    pack_send_buffer<parallelDim>(send_hv, send_dv, send_buffer, send_size, on_device);

    isend_and_wait_for_receive<DATA_TYPE>(tag, recv_counts_init, recv_req, send_displs, send_counts, send_req,
                                          send_buffer);

    /// Unpack
    unpack_recv_buffer<parallelDim>(recv_buffer, recv_size, recv_hv, recv_dv, on_device);

    wait_for_send(send_counts_init, send_req);

    deallocate_buffer<DATA_TYPE>(send_buffer, send_size, on_device);
    deallocate_buffer<DATA_TYPE>(recv_buffer, recv_size, on_device);

}

template <typename DATA_TYPE>
DATA_TYPE* allocate_buffer(const int buffer_size, const bool on_device) const {
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
void deallocate_buffer(DATA_TYPE* buffer, const int buffer_size, const bool on_device) const {
    if (on_device) {
        util::delete_devicemem(buffer, buffer_size);
    }
    else {
        util::delete_hostmem(buffer, buffer_size);
    }
}

void counts_displs_setup(const idx_t var_size, std::vector<int>& send_counts_init,
                         std::vector<int>& recv_counts_init, std::vector<int>& send_counts,
                         std::vector<int>& recv_counts, std::vector<int>& send_displs,
                         std::vector<int>& recv_displs) const {
    for (size_t jproc = 0; jproc < static_cast<size_t>(mpi_size_); ++jproc) {
        send_counts_init[jproc] = sendcounts_[jproc];
        recv_counts_init[jproc] = recvcounts_[jproc];
        send_counts[jproc]      = sendcounts_[jproc] * var_size;
        recv_counts[jproc]      = recvcounts_[jproc] * var_size;
        send_displs[jproc]      = senddispls_[jproc] * var_size;
        recv_displs[jproc]      = recvdispls_[jproc] * var_size;
    }
}
template <typename DATA_TYPE>
void ireceive(int tag, std::vector<int>& recv_displs, std::vector<int>& recv_counts,
              std::vector<eckit::mpi::Request>& recv_req, DATA_TYPE* recv_buffer) const {
    ATLAS_TRACE_MPI(IRECEIVE) {
        /// Let MPI know what we like to receive
        for (size_t jproc = 0; jproc < static_cast<size_t>(mpi_size_); ++jproc) {
            if (recv_counts[jproc] > 0) {
                recv_req[jproc] =
                    comm().iReceive(recv_buffer+recv_displs[jproc], recv_counts[jproc], jproc, tag);
            }
        }
    }
}

template <int ParallelDim, typename DATA_TYPE, int RANK>
void pack_send_buffer(ATLAS_MAYBE_UNUSED const array::ArrayView<DATA_TYPE, RANK>& hfield,
                                    const array::ArrayView<DATA_TYPE, RANK>& dfield, DATA_TYPE* send_buffer,
                                    int send_size, ATLAS_MAYBE_UNUSED const bool on_device) const {
    ATLAS_TRACE();
#if ATLAS_HAVE_GPU
    if (on_device) {
        halo_packer_hic<ParallelDim, DATA_TYPE, RANK>::pack(sendcnt_, sendmap_, hfield, dfield, send_buffer,
                                                             send_size);
    }
    else
#endif
        halo_packer<ParallelDim, RANK>::pack(sendcnt_, sendmap_, dfield, send_buffer, send_size);
}

template <int ParallelDim, typename DATA_TYPE, int RANK>
void unpack_recv_buffer(const DATA_TYPE* recv_buffer, int recv_size,
                                      ATLAS_MAYBE_UNUSED array::ArrayView<DATA_TYPE, RANK>& hfield,
                                      array::ArrayView<DATA_TYPE, RANK>& dfield,
                                      ATLAS_MAYBE_UNUSED const bool on_device) const {
    ATLAS_TRACE();
#if ATLAS_HAVE_GPU
    if (on_device) {
        halo_packer_hic<ParallelDim, DATA_TYPE, RANK>::unpack(recvcnt_, recvmap_, recv_buffer, recv_size, hfield,
                                                               dfield);
    }
    else
#endif
        halo_packer<ParallelDim, RANK>::unpack(recvcnt_, recvmap_, recv_buffer, recv_size, dfield);
}

template <typename DATA_TYPE>
void isend_and_wait_for_receive(int tag, std::vector<int>& recv_counts_init,
                                std::vector<eckit::mpi::Request>& recv_req, std::vector<int>& send_displs,
                                std::vector<int>& send_counts, std::vector<eckit::mpi::Request>& send_req,
                                DATA_TYPE* send_buffer) const {
    /// Send
    ATLAS_TRACE_MPI(ISEND) {
        for (size_t jproc = 0; jproc < static_cast<size_t>(mpi_size_); ++jproc) {
            if (send_counts[jproc] > 0) {
                send_req[jproc] = comm().iSend(send_buffer+send_displs[jproc], send_counts[jproc], jproc, tag);
            }
        }
    }

    /// Wait for receiving to finish
    ATLAS_TRACE_MPI(WAIT, "mpi-wait receive") {
        for (size_t jproc = 0; jproc < static_cast<size_t>(mpi_size_); ++jproc) {
            if (recv_counts_init[jproc] > 0) {
                comm().wait(recv_req[jproc]);
            }
        }
    }
}

void wait_for_send(std::vector<int>& send_counts_init, std::vector<eckit::mpi::Request>& send_req) const {
    ATLAS_TRACE_MPI(WAIT, "mpi-wait send") {
        for (size_t jproc = 0; jproc < static_cast<size_t>(mpi_size_); ++jproc) {
            if (send_counts_init[jproc] > 0) {
                comm().wait(send_req[jproc]);
            }
        }
    }
}

};

} // namespace atlas::parallel
