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
#include <iostream>
#include <numeric>
#include <sstream>
#include <stdexcept>

#include "eckit/log/Bytes.h"

#include "atlas/array.h"
#include "atlas/array/ArrayView.h"
#include "atlas/parallel/GatherScatter.h"
#include "atlas/parallel/mpi/Statistics.h"
#include "atlas/runtime/Log.h"
#include "atlas/runtime/Trace.h"

#include "atlas/parallel/omp/sort.h"

namespace atlas {
namespace parallel {

namespace {
struct IsGhostPoint {
    IsGhostPoint(const int mypart, const int part[], const idx_t ridx[], const idx_t base, const int N) {
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

struct Node {
    int p;
    idx_t i;
    gidx_t g;

    Node() = default;
    Node(gidx_t gid, int part, idx_t idx) {
        g = gid;
        p = part;
        i = idx;
    }

    bool operator<(const Node& other) const { return (g < other.g); }

    bool operator==(const Node& other) const { return (g == other.g); }
};

}  // namespace

GatherScatter::GatherScatter(): name_(), is_setup_(false) {
}

GatherScatter::GatherScatter(const std::string& name): name_(name), is_setup_(false) {
}

void GatherScatter::setup(const int part[], const idx_t remote_idx[], const int base, const gidx_t glb_idx[],
                          const int mask[], const idx_t parsize) {
    setup(mpi::comm().name(), part, remote_idx, base, glb_idx, parsize);
}

void GatherScatter::setup(const std::string& mpi_comm, const int part[], const idx_t remote_idx[], const int base, const gidx_t glb_idx[],
                          const int mask[], const idx_t parsize) {
    ATLAS_TRACE("GatherScatter::setup");
    comm_ = &mpi::comm(mpi_comm);
    myproc = comm().rank();
    nproc  = comm().size();
    parsize_ = parsize;

    glbcounts_.resize(nproc);
    glbcounts_.assign(nproc, 0);
    glbdispls_.resize(nproc);
    glbdispls_.assign(nproc, 0);

    std::vector<gidx_t> sendnodes_gidx(parsize_);
    std::vector<int> sendnodes_part(parsize_);
    std::vector<idx_t> sendnodes_ridx(parsize_);

    loccnt_ = 0;
    for (idx_t n = 0; n < parsize_; ++n) {
        if (!mask[n]) {
            sendnodes_gidx[loccnt_] = glb_idx[n];
            sendnodes_part[loccnt_] = part[n];
            sendnodes_ridx[loccnt_] = remote_idx[n] - base;
            ++loccnt_;
        }
    }

    ATLAS_TRACE_MPI(ALLGATHER) {
        comm().allGather(loccnt_, glbcounts_.begin(), glbcounts_.end());
    }

    size_t glbcnt_size_t = std::accumulate(glbcounts_.begin(), glbcounts_.end(), size_t(0));
    if (glbcnt_size_t > std::numeric_limits<int>::max()) {
        ATLAS_THROW_EXCEPTION("Due to limitation of MPI we cannot use larger counts");
    }

    glbcnt_ = std::accumulate(glbcounts_.begin(), glbcounts_.end(), size_t(0));

    glbdispls_[0] = 0;
    for (idx_t jproc = 1; jproc < nproc; ++jproc)  // start at 1
    {
        glbdispls_[jproc] = glbcounts_[jproc - 1] + glbdispls_[jproc - 1];
    }

    idx_t nb_recv_nodes = glbcnt_;
    std::vector<Node> node_sort;

    try {
        node_sort.resize(nb_recv_nodes);
    }
    catch (const std::bad_alloc& e) {
        ATLAS_THROW_EXCEPTION("Could not allocate node_sort with size " << eckit::Bytes(nb_recv_nodes * sizeof(Node)));
    }

    {
        std::vector<gidx_t> recvnodes_gidx;
        try {
            recvnodes_gidx.resize(glbcnt_);
        }
        catch (const std::bad_alloc& e) {
            ATLAS_THROW_EXCEPTION("Could not allocate recvnodes_gidx with size "
                                  << eckit::Bytes(glbcnt_ * sizeof(gidx_t)));
        }
        ATLAS_TRACE_MPI(ALLGATHER) {
            comm().allGatherv(sendnodes_gidx.begin(), sendnodes_gidx.begin() + loccnt_, recvnodes_gidx.data(),
                              glbcounts_.data(), glbdispls_.data());
        }
        atlas_omp_parallel_for(idx_t n = 0; n < nb_recv_nodes; ++n) {
            node_sort[n].g = recvnodes_gidx[n];
        }
        sendnodes_gidx.clear();
    }


    {
        std::vector<int> recvnodes_part;
        try {
            recvnodes_part.resize(glbcnt_);
        }
        catch (const std::bad_alloc& e) {
            ATLAS_THROW_EXCEPTION("Could not allocate recvnodes_part with size "
                                  << eckit::Bytes(glbcnt_ * sizeof(int)));
        }
        ATLAS_TRACE_MPI(ALLGATHER) {
            comm().allGatherv(sendnodes_part.begin(), sendnodes_part.begin() + loccnt_, recvnodes_part.data(),
                              glbcounts_.data(), glbdispls_.data());
        }
        atlas_omp_parallel_for(idx_t n = 0; n < nb_recv_nodes; ++n) {
            node_sort[n].p = recvnodes_part[n];
        }
        sendnodes_part.clear();
    }

    {
        std::vector<idx_t> recvnodes_ridx;
        try {
            recvnodes_ridx.resize(glbcnt_);
        }
        catch (const std::bad_alloc& e) {
            ATLAS_THROW_EXCEPTION("Could not allocate recvnodes_ridx with size "
                                  << eckit::Bytes(glbcnt_ * sizeof(idx_t)));
        }
        ATLAS_TRACE_MPI(ALLGATHER) {
            comm().allGatherv(sendnodes_ridx.begin(), sendnodes_ridx.begin() + loccnt_, recvnodes_ridx.data(),
                              glbcounts_.data(), glbdispls_.data());
        }
        atlas_omp_parallel_for(idx_t n = 0; n < nb_recv_nodes; ++n) {
            node_sort[n].i = recvnodes_ridx[n];
        }
        sendnodes_ridx.clear();
    }


    // Sort on "g" member, and remove duplicates
    ATLAS_TRACE_SCOPE("sorting") {
        //        omp::sort(node_sort.begin(), node_sort.end());
        std::sort(node_sort.begin(), node_sort.end());
        node_sort.erase(std::unique(node_sort.begin(), node_sort.end()), node_sort.end());
    }

    glbcounts_.assign(nproc, 0);
    glbdispls_.assign(nproc, 0);

    for (size_t n = 0; n < node_sort.size(); ++n) {
        ++glbcounts_[node_sort[n].p];
    }

    glbdispls_[0] = 0;

    for (idx_t jproc = 1; jproc < nproc; ++jproc)  // start at 1
    {
        glbdispls_[jproc] = glbcounts_[jproc - 1] + glbdispls_[jproc - 1];
    }

    glbcnt_ = std::accumulate(glbcounts_.begin(), glbcounts_.end(), size_t(0));
    loccnt_ = glbcounts_[myproc];

    glbmap_.clear();
    glbmap_.resize(glbcnt_);
    locmap_.clear();
    locmap_.resize(loccnt_);
    std::vector<int> idx(nproc, 0);

    size_t n{0};
    for (const auto& node : node_sort) {
        idx_t jproc                             = node.p;
        glbmap_[glbdispls_[jproc] + idx[jproc]] = n++;

        if (jproc == myproc) {
            locmap_[idx[jproc]] = node.i;
        }

        ++idx[jproc];
    }

    is_setup_ = true;
}

void GatherScatter::setup(const int part[], const idx_t remote_idx[], const int base, const gidx_t glb_idx[],
                          const idx_t parsize) {
    setup(mpi::comm().name(), part, remote_idx, base, glb_idx, parsize);
}

void GatherScatter::setup(const std::string& mpi_comm, const int part[], const idx_t remote_idx[], const int base, const gidx_t glb_idx[],
                          const idx_t parsize) {
    std::vector<int> mask(parsize);
    {
        int mypart = mpi::comm(mpi_comm).rank();
        IsGhostPoint is_ghost(mypart, part, remote_idx, base, parsize);
        for (idx_t jj = 0; jj < parsize; ++jj) {
            mask[jj] = is_ghost(jj) ? 1 : 0;
        }
    }
    setup(mpi_comm, part, remote_idx, base, glb_idx, mask.data(), parsize);
}

/////////////////////

GatherScatter* atlas__GatherScatter__new() {
    return new GatherScatter();
}

void atlas__GatherScatter__delete(GatherScatter* This) {
    delete This;
}

void atlas__GatherScatter__setup32(GatherScatter* This, int part[], idx_t remote_idx[], int base, int glb_idx[],
                                   int parsize) {
#if ATLAS_BITS_GLOBAL == 32
    This->setup(part, remote_idx, base, glb_idx, parsize);
#else
    std::vector<gidx_t> glb_idx_convert(parsize);
    for (int j = 0; j < parsize; ++j) {
        glb_idx_convert[j] = glb_idx[j];
    }
    This->setup(part, remote_idx, base, glb_idx_convert.data(), parsize);
#endif
}

void atlas__GatherScatter__setup64(GatherScatter* This, int part[], idx_t remote_idx[], int base, long glb_idx[],
                                   int parsize) {
#if ATLAS_BITS_GLOBAL == 64
    This->setup(part, remote_idx, base, glb_idx, parsize);
#else
    std::vector<gidx_t> glb_idx_convert(parsize);
    for (idx_t j = 0; j < parsize; ++j) {
        glb_idx_convert[j] = glb_idx[j];
    }
    This->setup(part, remote_idx, base, glb_idx_convert.data(), parsize);
#endif
}

int atlas__GatherScatter__glb_dof(GatherScatter* This) {
    return This->glb_dof();
}

void atlas__GatherScatter__gather_int(GatherScatter* This, int lfield[], int lvar_strides[], int lvar_extents[],
                                      int lvar_rank, int gfield[], int gvar_strides[], int gvar_extents[],
                                      int gvar_rank) {
    std::vector<idx_t> lvstrides(lvar_rank);
    std::vector<idx_t> lvextents(lvar_rank);
    std::vector<idx_t> gvstrides(gvar_rank);
    std::vector<idx_t> gvextents(gvar_rank);
    for (int n = 0; n < lvar_rank; ++n) {
        lvstrides[n] = lvar_strides[n];
        lvextents[n] = lvar_extents[n];
    }
    for (int n = 0; n < gvar_rank; ++n) {
        gvstrides[n] = gvar_strides[n];
        gvextents[n] = gvar_extents[n];
    }
    This->gather(lfield, lvstrides.data(), lvextents.data(), lvar_rank, gfield, gvstrides.data(), gvextents.data(),
                 gvar_rank);
}

void atlas__GatherScatter__gather_long(GatherScatter* This, long lfield[], int lvar_strides[], int lvar_extents[],
                                       int lvar_rank, long gfield[], int gvar_strides[], int gvar_extents[],
                                       int gvar_rank) {
    std::vector<idx_t> lvstrides(lvar_rank);
    std::vector<idx_t> lvextents(lvar_rank);
    std::vector<idx_t> gvstrides(gvar_rank);
    std::vector<idx_t> gvextents(gvar_rank);
    for (int n = 0; n < lvar_rank; ++n) {
        lvstrides[n] = lvar_strides[n];
        lvextents[n] = lvar_extents[n];
    }
    for (int n = 0; n < gvar_rank; ++n) {
        gvstrides[n] = gvar_strides[n];
        gvextents[n] = gvar_extents[n];
    }
    This->gather(lfield, lvstrides.data(), lvextents.data(), lvar_rank, gfield, gvstrides.data(), gvextents.data(),
                 gvar_rank);
}

void atlas__GatherScatter__gather_float(GatherScatter* This, float lfield[], int lvar_strides[], int lvar_extents[],
                                        int lvar_rank, float gfield[], int gvar_strides[], int gvar_extents[],
                                        int gvar_rank) {
    std::vector<idx_t> lvstrides(lvar_rank);
    std::vector<idx_t> lvextents(lvar_rank);
    std::vector<idx_t> gvstrides(gvar_rank);
    std::vector<idx_t> gvextents(gvar_rank);
    for (int n = 0; n < lvar_rank; ++n) {
        lvstrides[n] = lvar_strides[n];
        lvextents[n] = lvar_extents[n];
    }
    for (int n = 0; n < gvar_rank; ++n) {
        gvstrides[n] = gvar_strides[n];
        gvextents[n] = gvar_extents[n];
    }
    This->gather(lfield, lvstrides.data(), lvextents.data(), lvar_rank, gfield, gvstrides.data(), gvextents.data(),
                 gvar_rank);
}

void atlas__GatherScatter__gather_double(GatherScatter* This, double lfield[], int lvar_strides[], int lvar_extents[],
                                         int lvar_rank, double gfield[], int gvar_strides[], int gvar_extents[],
                                         int gvar_rank) {
    std::vector<idx_t> lvstrides(lvar_rank);
    std::vector<idx_t> lvextents(lvar_rank);
    std::vector<idx_t> gvstrides(gvar_rank);
    std::vector<idx_t> gvextents(gvar_rank);
    for (int n = 0; n < lvar_rank; ++n) {
        lvstrides[n] = lvar_strides[n];
        lvextents[n] = lvar_extents[n];
    }
    for (int n = 0; n < gvar_rank; ++n) {
        gvstrides[n] = gvar_strides[n];
        gvextents[n] = gvar_extents[n];
    }
    This->gather(lfield, lvstrides.data(), lvextents.data(), lvar_rank, gfield, gvstrides.data(), gvextents.data(),
                 gvar_rank);
}

void atlas__GatherScatter__scatter_int(GatherScatter* This, int gfield[], int gvar_strides[], int gvar_extents[],
                                       int gvar_rank, int lfield[], int lvar_strides[], int lvar_extents[],
                                       int lvar_rank) {
    std::vector<idx_t> lvstrides(lvar_rank);
    std::vector<idx_t> lvextents(lvar_rank);
    std::vector<idx_t> gvstrides(gvar_rank);
    std::vector<idx_t> gvextents(gvar_rank);
    for (int n = 0; n < lvar_rank; ++n) {
        lvstrides[n] = lvar_strides[n];
        lvextents[n] = lvar_extents[n];
    }
    for (int n = 0; n < gvar_rank; ++n) {
        gvstrides[n] = gvar_strides[n];
        gvextents[n] = gvar_extents[n];
    }
    This->scatter(gfield, gvstrides.data(), gvextents.data(), gvar_rank, lfield, lvstrides.data(), lvextents.data(),
                  lvar_rank);
}

void atlas__GatherScatter__scatter_long(GatherScatter* This, long gfield[], int gvar_strides[], int gvar_extents[],
                                        int gvar_rank, long lfield[], int lvar_strides[], int lvar_extents[],
                                        int lvar_rank) {
    std::vector<idx_t> lvstrides(lvar_rank);
    std::vector<idx_t> lvextents(lvar_rank);
    std::vector<idx_t> gvstrides(gvar_rank);
    std::vector<idx_t> gvextents(gvar_rank);
    for (int n = 0; n < lvar_rank; ++n) {
        lvstrides[n] = lvar_strides[n];
        lvextents[n] = lvar_extents[n];
    }
    for (int n = 0; n < gvar_rank; ++n) {
        gvstrides[n] = gvar_strides[n];
        gvextents[n] = gvar_extents[n];
    }
    This->scatter(gfield, gvstrides.data(), gvextents.data(), gvar_rank, lfield, lvstrides.data(), lvextents.data(),
                  lvar_rank);
}

void atlas__GatherScatter__scatter_float(GatherScatter* This, float gfield[], int gvar_strides[], int gvar_extents[],
                                         int gvar_rank, float lfield[], int lvar_strides[], int lvar_extents[],
                                         int lvar_rank) {
    std::vector<idx_t> lvstrides(lvar_rank);
    std::vector<idx_t> lvextents(lvar_rank);
    std::vector<idx_t> gvstrides(gvar_rank);
    std::vector<idx_t> gvextents(gvar_rank);
    for (int n = 0; n < lvar_rank; ++n) {
        lvstrides[n] = lvar_strides[n];
        lvextents[n] = lvar_extents[n];
    }
    for (int n = 0; n < gvar_rank; ++n) {
        gvstrides[n] = gvar_strides[n];
        gvextents[n] = gvar_extents[n];
    }
    This->scatter(gfield, gvstrides.data(), gvextents.data(), gvar_rank, lfield, lvstrides.data(), lvextents.data(),
                  lvar_rank);
}

void atlas__GatherScatter__scatter_double(GatherScatter* This, double gfield[], int gvar_strides[], int gvar_extents[],
                                          int gvar_rank, double lfield[], int lvar_strides[], int lvar_extents[],
                                          int lvar_rank) {
    std::vector<idx_t> lvstrides(lvar_rank);
    std::vector<idx_t> lvextents(lvar_rank);
    std::vector<idx_t> gvstrides(gvar_rank);
    std::vector<idx_t> gvextents(gvar_rank);
    for (int n = 0; n < lvar_rank; ++n) {
        lvstrides[n] = lvar_strides[n];
        lvextents[n] = lvar_extents[n];
    }
    for (int n = 0; n < gvar_rank; ++n) {
        gvstrides[n] = gvar_strides[n];
        gvextents[n] = gvar_extents[n];
    }
    This->scatter(gfield, gvstrides.data(), gvextents.data(), gvar_rank, lfield, lvstrides.data(), lvextents.data(),
                  lvar_rank);
}

/////////////////////

}  // namespace parallel
}  // namespace atlas
