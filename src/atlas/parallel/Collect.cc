/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/parallel/Collect.h"


namespace atlas::parallel {

void Collect::setup(const std::string& mpi_comm,
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

    atlas_omp_parallel {
        std::vector<int> recvcounts(nproc, 0);
        atlas_omp_for(std::size_t jrecv = 0; jrecv < recvcnt_; ++jrecv) {
            int p = recv_part[jrecv];
            ++recvcounts[p];
        }
        atlas_omp_critical {
            for (int p=0; p < nproc; ++p) {
                recvcounts_[p] += recvcounts[p];
            }
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
    Log::trace() << __LINE__ << std::endl;

    for (int jrecv = 0; jrecv < recvcnt_; ++jrecv) {
        const int p            = recv_part[jrecv];
        const int req_idx      = recvdispls_[p] + cnt[p];
        send_requests[req_idx] = recv_remote_idx[jrecv] - recv_idx_base;
        recvmap_[req_idx]      = jrecv;
        cnt[p]++;
    }

    Log::trace() << __LINE__ << std::endl;

    /*
    Fill vector "recv_requests" with what is needed by other procs
    */
    ATLAS_TRACE_MPI(ALLTOALL) {
        ATLAS_ASSERT(send_requests.size() == recvcnt_);
        ATLAS_ASSERT(std::accumulate(recvcounts_.begin(),recvcounts_.end(),0) == recvcnt_);
        ATLAS_ASSERT(recv_requests.size() == sendcnt_);
        ATLAS_ASSERT(std::accumulate(sendcounts_.begin(),sendcounts_.end(),0) == sendcnt_);

        comm().allToAllv(send_requests.data(), recvcounts_.data(), recvdispls_.data(), recv_requests.data(),
                              sendcounts_.data(), senddispls_.data());
    }
    Log::trace() << __LINE__ << std::endl;


    /*
    What needs to be sent to other procs is asked by remote_idx, which is local
    here
    */
    sendmap_.resize(sendcnt_);
    for (int jj = 0; jj < sendcnt_; ++jj) {
        sendmap_[jj] = recv_requests[jj];
    }
    Log::trace() << __LINE__ << std::endl;

    is_setup_        = true;
}


void Collect::counts_displs_setup(const idx_t var_size, std::vector<int>& send_counts_init,
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

void Collect::send_counts_displs_setup(const idx_t var_size, std::vector<int>& send_counts_init,
                            std::vector<int>& send_counts,
                            std::vector<int>& send_displs) const {
    for (size_t jproc = 0; jproc < static_cast<size_t>(mpi_size_); ++jproc) {
        send_counts_init[jproc] = sendcounts_[jproc];
        send_counts[jproc]      = sendcounts_[jproc] * var_size;
        send_displs[jproc]      = senddispls_[jproc] * var_size;
    }
}



void Collect::wait_for_send(std::vector<int>& send_counts_init, std::vector<eckit::mpi::Request>& send_req) const {
    ATLAS_TRACE_MPI(WAIT, "mpi-wait send") {
        for (size_t jproc = 0; jproc < static_cast<size_t>(mpi_size_); ++jproc) {
            if (send_counts_init[jproc] > 0) {
                comm().wait(send_req[jproc]);
            }
        }
    }
}

} // namespace atlas::parallel
