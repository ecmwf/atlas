/*
 * (C) Copyright 2025 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "Locate.h"

#include <unordered_map>

#include "atlas/parallel/mpi/mpi.h"
#include "atlas/runtime/Trace.h"

namespace atlas::parallel {

// Given distribution and global_index, find the corresponding partition for each global index
// The global-index is typically 1-based
void Locator::locate_partition(
    fspan<const int> distribution, const int distribution_base,
    span<const gidx_t> global_index, const gidx_t global_index_base,
    span<int> partition, const int partition_base) {

    ATLAS_TRACE("atlas::util::locate_partition");
    std::size_t size = global_index.size();
    for (std::size_t j = 0; j < size; ++j) {
        auto gidx = global_index[j] - global_index_base;
        ATLAS_ASSERT(gidx < distribution.size());
        auto p = distribution[gidx] - distribution_base;
        partition[j] = p + partition_base;
    }
}

// Given distribution and global_index, find the corresponding partition for each global index
// The global-index is typically 1-based
void Locator::locate_partition(
    span<const int> distribution, const int distribution_base,
    span<const gidx_t> global_index, const gidx_t global_index_base,
    span<int> partition, const int partition_base) {
    ATLAS_TRACE("atlas::util::locate_partition");
    std::size_t size = global_index.size();
    for (std::size_t j = 0; j < size; ++j) {
        gidx_t gidx = global_index[j] - global_index_base;
        ATLAS_ASSERT(gidx < distribution.size());
        int p = distribution[gidx] - distribution_base;
        partition[j] = p + partition_base;
    }
}

// Given global_index, find the corresponding partition and remote_index
// This could be costly as it involves memory and communication
void Locator::locate_remote_index(
    const std::string_view mpi_comm,
    span<const gidx_t> my_glb_idx, const gidx_t my_glb_idx_base, span<const int> my_ghost,
    span<const gidx_t> global_index, const gidx_t global_index_base,
    span<const int> partition, const int partition_base,
    span<idx_t> remote_index, const idx_t remote_index_base) {

    ATLAS_TRACE("atlas::util::locate_remote_index");
    auto& comm = atlas::mpi::comm(mpi_comm);
    int mpi_size = comm.size();
    std::size_t my_size = my_glb_idx.size();
    std::size_t size = global_index.size();
    std::vector<std::vector<gidx_t>> recv_gidx(mpi_size);
    {
        std::vector<std::vector<gidx_t>> send_gidx(mpi_size);
        std::vector<std::size_t> send_counts(mpi_size,0);
        for (std::size_t j=0; j<size; ++j) {
            int p = partition[j] - partition_base;
            ++send_counts[p];
        }
        for (int p=0; p<mpi_size; ++p) {
            send_gidx[p].reserve(send_counts[p]);
        }
        for (std::size_t j=0; j<size; ++j) {
            gidx_t gidx = global_index[j] - global_index_base;
            int p = partition[j] - partition_base;
            send_gidx[p].emplace_back(gidx);
        }
        comm.allToAll(send_gidx, recv_gidx);
    }

    std::vector<std::vector<idx_t>> recv_ridx(mpi_size);
    {
        std::vector<std::vector<idx_t>> send_ridx(mpi_size);

        std::unordered_map<gidx_t, idx_t> gidx_to_ridx;
        for (idx_t i=0; i<my_size; ++i) {
            if (not my_ghost[i]) {
                gidx_to_ridx[my_glb_idx[i] - my_glb_idx_base] = i;
            }
        }
        for( int p=0; p < recv_gidx.size(); ++p ) {
            for( auto gidx : recv_gidx[p] ) {
                auto it = gidx_to_ridx.find(gidx);
                ATLAS_ASSERT( it != gidx_to_ridx.end() );
                send_ridx[p].emplace_back(it->second);
            }
        }
        comm.allToAll(send_ridx, recv_ridx);
    }
    std::vector<std::size_t> c(mpi_size,0);
    for( idx_t j=0; j<size; ++j) {
        int p = partition[j] - partition_base;
        remote_index[j] = recv_ridx[p][c[p]++] + remote_index_base;
    }
}


// Given global_index, find the corresponding partition and remote_index
// This could be costly as it involves memory and communication
// local information: my_size, my_glb_idx, my_ghost
// global information: distribution
// requested indices to locate: size, global_index
// output of locate: partition, remote_index, remote_index_base
void Locator::locate( const std::string_view mpi_comm, 
    span<const gidx_t> my_glb_idx, const gidx_t my_glb_idx_base, span<const int> my_ghost,
    fspan<const int> distribution, const int distribution_base,
    span<const gidx_t> global_index, const gidx_t global_index_base,
    span<int> partition, const int partition_base,
    span<idx_t> remote_index, const idx_t remote_index_base) {

    ATLAS_TRACE("atlas::util::locate");

    locate_partition(distribution, distribution_base, global_index, global_index_base,
                     partition, partition_base);

    locate_remote_index(mpi_comm, my_glb_idx, my_glb_idx_base, my_ghost,
                        global_index, global_index_base,
                        partition, partition_base,
                        remote_index, remote_index_base);
}

void Locator::locate( const std::string_view mpi_comm, 
    span<const gidx_t> my_glb_idx, const gidx_t my_glb_idx_base, span<const int> my_ghost,
    span<const int> distribution, const int distribution_base,
    span<const gidx_t> global_index, const gidx_t global_index_base,
    span<int> partition, const int partition_base,
    span<idx_t> remote_index, const idx_t remote_index_base) {

    ATLAS_TRACE("atlas::util::locate");

    locate_partition(distribution, distribution_base,
                     global_index, global_index_base,
                     partition, partition_base);

    locate_remote_index(mpi_comm, my_glb_idx, my_glb_idx_base, my_ghost,
                        global_index, global_index_base,
                        partition, partition_base,
                        remote_index, remote_index_base);
}

} // atlas::parallel
