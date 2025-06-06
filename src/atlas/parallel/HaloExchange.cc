/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

/// @file HaloExchange.cc
/// @author Willem Deconinck
/// @date   Nov 2013

#include <memory>
#include <numeric>
#include <sstream>
#include <stdexcept>

#include "pluto/pluto.h"

#include "atlas/array/Array.h"
#include "atlas/parallel/HaloExchange.h"
#include "atlas/parallel/mpi/Statistics.h"
#include "atlas/util/vector.h"

namespace atlas {
namespace parallel {

namespace {
struct IsGhostPoint {
    IsGhostPoint(const int mypart, const int part[], const idx_t ridx[], const idx_t base, const int /*N*/) {
        part_   = part;
        ridx_   = ridx;
        base_   = base;
        mypart_ = mypart;
    }

    bool operator()(idx_t idx) {
        if (part_[idx] != mypart_) {
            return true;
        }
        if (ridx_[idx] != base_ + idx) {
            return true;
        }
        return false;
    }
    int mypart_;
    const int* part_;
    const idx_t* ridx_;
    idx_t base_;
};
}  // namespace

HaloExchange::HaloExchange() :
    HaloExchange("") {
}

HaloExchange::HaloExchange(const std::string& name):
    name_(name),
    is_setup_(false) {
}

HaloExchange::~HaloExchange() = default;

void HaloExchange::setup(const int part[], const idx_t remote_idx[], const int base, const idx_t size) {
    setup(mpi::comm().name(), part, remote_idx, base, size);
}

void HaloExchange::setup(const int part[], const idx_t remote_idx[], const int base, idx_t parsize, idx_t halo_begin) {
    setup(mpi::comm().name(), part, remote_idx, base, parsize, halo_begin);
}

void HaloExchange::setup(const std::string& mpi_comm, const int part[], const idx_t remote_idx[], const int base, const idx_t size) {
    setup(mpi_comm, part, remote_idx, base, size, 0);
}

void HaloExchange::setup(const std::string& mpi_comm, const int part[], const idx_t remote_idx[], const int base, idx_t parsize, idx_t halo_begin) {
    ATLAS_TRACE("HaloExchange::setup");
    comm_ = &mpi::comm(mpi_comm);
    myproc = comm().rank();
    nproc  = comm().size();

    parsize_ = parsize;
    sendcounts_.resize(nproc);
    sendcounts_.assign(nproc, 0);
    recvcounts_.resize(nproc);
    recvcounts_.assign(nproc, 0);
    senddispls_.resize(nproc);
    senddispls_.assign(nproc, 0);
    recvdispls_.resize(nproc);
    recvdispls_.assign(nproc, 0);

    /*
    Find the amount of nodes this proc has to receive from each other proc
    */

    IsGhostPoint is_ghost(myproc, part, remote_idx, base, parsize_);
    atlas::vector<idx_t> ghost_points(parsize_);
    idx_t nghost = 0;

    atlas_omp_parallel_for(int jj = halo_begin; jj < parsize_; ++jj) {
        if (is_ghost(jj)) {
            int p = part[jj];
            atlas_omp_critical {
                ++recvcounts_[p];
                ghost_points[nghost] = jj;
                nghost++;
            }
        }
    }

    recvcnt_ = std::accumulate(recvcounts_.begin(), recvcounts_.end(), 0);

    /*
    Find the amount of nodes this proc has to send to each other proc
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
    for (idx_t jghost = 0; jghost < nghost; ++jghost) {
        const idx_t jj         = ghost_points[jghost];
        const int p            = part[jj];
        const int req_idx      = recvdispls_[p] + cnt[p];
        send_requests[req_idx] = remote_idx[jj] - base;
        recvmap_[req_idx]      = jj;
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
    backdoor.parsize = parsize_;
}

void HaloExchange::wait_for_send(std::vector<int>& send_counts_init, std::vector<eckit::mpi::Request>& send_req) const {
    ATLAS_TRACE_MPI(WAIT, "mpi-wait send") {
        for (int jproc = 0; jproc < nproc; ++jproc) {
            if (send_counts_init[jproc] > 0) {
                comm().wait(send_req[jproc]);
            }
        }
    }
}

int HaloExchange::devices() const {
    return pluto::devices();
}

namespace {

template <typename Value>
void execute_halo_exchange(HaloExchange* This, Value field[], int var_strides[], int var_extents[], int var_rank) {
    // WARNING: Only works if there is only one parallel dimension AND being
    // slowest moving

    array::ArrayShape shape{This->backdoor.parsize};
    for (int j = 0; j < var_rank; ++j) {
        shape.push_back(var_extents[j]);
    }

    array::ArrayStrides strides{var_extents[0] * var_strides[0]};
    for (int j = 0; j < var_rank; ++j) {
        strides.push_back(var_strides[j]);
    }

    std::unique_ptr<array::Array> arr(array::Array::wrap(field, array::ArraySpec{shape, strides}));

    switch (arr->rank()) {
        case 1: {
            This->execute<Value, 1>(*arr);
            break;
        }
        case 2: {
            This->execute<Value, 2>(*arr);
            break;
        }
        case 3: {
            This->execute<Value, 3>(*arr);
            break;
        }
        case 4: {
            This->execute<Value, 4>(*arr);
            break;
        }
        default:
            throw_NotImplemented("Rank not supported in halo exchange", Here());
    }
}

template <typename Value>
void execute_adjoint_halo_exchange(HaloExchange* This, Value field[], int var_strides[], int var_extents[],
                                   int var_rank) {
    // WARNING: Only works if there is only one parallel dimension AND being
    // slowest moving

    array::ArrayShape shape{This->backdoor.parsize};
    for (int j = 0; j < var_rank; ++j) {
        shape.push_back(var_extents[j]);
    }

    array::ArrayStrides strides{var_extents[0] * var_strides[0]};
    for (int j = 0; j < var_rank; ++j) {
        strides.push_back(var_strides[j]);
    }

    std::unique_ptr<array::Array> arr(array::Array::wrap(field, array::ArraySpec{shape, strides}));

    switch (arr->rank()) {
        case 1: {
            This->execute_adjoint<Value, 1>(*arr);
            break;
        }
        case 2: {
            This->execute_adjoint<Value, 2>(*arr);
            break;
        }
        case 3: {
            This->execute_adjoint<Value, 3>(*arr);
            break;
        }
        case 4: {
            This->execute_adjoint<Value, 4>(*arr);
            break;
        }
        default:
            throw_NotImplemented("Rank not supported in halo exchange", Here());
    }
}
}  // namespace

extern "C" {

HaloExchange* atlas__HaloExchange__new() {
    return new HaloExchange();
}

void atlas__HaloExchange__delete(HaloExchange* This) {
    delete This;
}

void atlas__HaloExchange__setup(HaloExchange* This, int part[], idx_t remote_idx[], int base, int size) {
    This->setup(part, remote_idx, base, size);
}


void atlas__HaloExchange__execute_adjoint_strided_int(HaloExchange* This, int field[], int var_strides[],
                                                      int var_extents[], int var_rank) {
    execute_adjoint_halo_exchange(This, field, var_strides, var_extents, var_rank);
}

void atlas__HaloExchange__execute_adjoint_strided_long(HaloExchange* This, long field[], int var_strides[],
                                                       int var_extents[], int var_rank) {
    execute_adjoint_halo_exchange(This, field, var_strides, var_extents, var_rank);
}

void atlas__HaloExchange__execute_adjoint_strided_float(HaloExchange* This, float field[], int var_strides[],
                                                        int var_extents[], int var_rank) {
    execute_adjoint_halo_exchange(This, field, var_strides, var_extents, var_rank);
}

void atlas__HaloExchange__execute_adjoint_strided_double(HaloExchange* This, double field[], int var_strides[],
                                                         int var_extents[], int var_rank) {
    execute_adjoint_halo_exchange(This, field, var_strides, var_extents, var_rank);
}

void atlas__HaloExchange__execute_strided_int(HaloExchange* This, int field[], int var_strides[], int var_extents[],
                                              int var_rank) {
    execute_halo_exchange(This, field, var_strides, var_extents, var_rank);
}

void atlas__HaloExchange__execute_strided_long(HaloExchange* This, long field[], int var_strides[], int var_extents[],
                                               int var_rank) {
    execute_halo_exchange(This, field, var_strides, var_extents, var_rank);
}

void atlas__HaloExchange__execute_strided_float(HaloExchange* This, float field[], int var_strides[], int var_extents[],
                                                int var_rank) {
    execute_halo_exchange(This, field, var_strides, var_extents, var_rank);
}

void atlas__HaloExchange__execute_strided_double(HaloExchange* This, double field[], int var_strides[],
                                                 int var_extents[], int var_rank) {
    execute_halo_exchange(This, field, var_strides, var_extents, var_rank);
}
}

/////////////////////

}  // namespace parallel
}  // namespace atlas
