/*
 * (C) Copyright 2025 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#pragma once

#include <vector>
#include <unordered_map>

#include "atlas/grid.h"
#include "atlas/functionspace.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/util/vector.h"

namespace atlas::util {


inline void locate_partition(const grid::Distribution& distribution, std::size_t size,
    const gidx_t global_index[], int partition[]) {

    ATLAS_TRACE("atlas::functionspace::locate_partition");
    for (std::size_t j = 0; j < size; ++j) {
        gidx_t gidx = global_index[j];
        ATLAS_ASSERT(gidx-1 < distribution.size());
        // assert_in_bounds(gidx-1, partition_global);
        int p = distribution.partition(gidx-1);
        partition[j] = p;
    }
}


// Given global_index, find the corresponding partition and remote_index
// This could be costly as it involves memory and communication
inline void locate_remote_index(
    const mpi::Comm& comm,
    std::size_t my_size, const gidx_t my_glb_idx[], const int my_ghost[],
    std::size_t size, const gidx_t global_index[], const int partition[],
    idx_t remote_index[], idx_t remote_index_base) {

    ATLAS_TRACE("atlas::functionspace::locate_remote_index");
    int mpi_size = comm.size();

    std::vector<std::vector<gidx_t>> recv_gidx(mpi_size);
    {
        std::vector<std::vector<gidx_t>> send_gidx(mpi_size);
        std::vector<std::size_t> send_counts(mpi_size,0);
        for (std::size_t j=0; j<size; ++j) {
            int p = partition[j];
            ++send_counts[p];
        }
        for (int p=0; p<mpi_size; ++p) {
            send_gidx[p].reserve(send_counts[p]);
        }
        for (std::size_t j=0; j<size; ++j) {
            gidx_t gidx = global_index[j];
            int p = partition[j];
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
                gidx_to_ridx[my_glb_idx[i]] = i;
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
        int p = partition[j];
        remote_index[j] = recv_ridx[p][c[p]++] + remote_index_base;
    }
}


// Given global_index, find the corresponding partition and remote_index
// This could be costly as it involves memory and communication
// local information: my_size, my_glb_idx, my_ghost
// global information: distribution
// requested indices to locate: size, global_index
// output of locate: partition, remote_index, remote_index_base
inline void locate( const mpi::Comm& comm, std::size_t my_size, const gidx_t my_glb_idx[], int my_ghost[],
    const grid::Distribution& distribution, std::size_t size, const gidx_t global_index[],
    int partition[], idx_t remote_index[], idx_t remote_index_base) {

    ATLAS_TRACE("atlas::functionspace::locate");
    int mpi_size = comm.size();

    std::vector<std::vector<gidx_t>> recv_gidx(mpi_size);
    {
        std::vector<std::vector<gidx_t>> send_gidx(mpi_size);
        std::vector<std::size_t> send_counts(mpi_size,0);
        for (std::size_t j=0; j<size; ++j) {
            gidx_t gidx = global_index[j];
            ATLAS_ASSERT(gidx-1 < distribution.size());
            // assert_in_bounds(gidx-1, partition_global);
            int p = distribution.partition(gidx-1);
            partition[j] = p;
            ++send_counts[p];
        }
        for (int p=0; p<mpi_size; ++p) {
            send_gidx[p].reserve(send_counts[p]);
        }
        for (std::size_t j=0; j<size; ++j) {
            gidx_t gidx = global_index[j];
            int p = partition[j];
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
                gidx_to_ridx[my_glb_idx[i]] = i;
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
        int p = partition[j];
        remote_index[j] = recv_ridx[p][c[p]++] + remote_index_base;
    }
}


inline void locate(const FunctionSpace& fs, const grid::Distribution& distribution, std::size_t size,
    const gidx_t global_index[], int partition[], idx_t remote_index[], idx_t remote_index_base) {
    ATLAS_TRACE("atlas::functionspace::locate");
    auto& comm = atlas::mpi::comm(fs.mpi_comm());
    auto glb_idx = array::make_view<gidx_t,1>(fs.global_index());
    auto ghost = array::make_view<int,1>(fs.ghost());
    locate(comm, glb_idx.size(), glb_idx.data(), ghost.data(), distribution, size,
        global_index, partition, remote_index, remote_index_base);
}


inline void locate(const FunctionSpace& fs, std::size_t size, const gidx_t global_index[], int partition[], idx_t remote_index[], idx_t remote_index_base) {
    ATLAS_TRACE("atlas::functionspace::locate");
    auto& comm = atlas::mpi::comm(fs.mpi_comm());
    int mpi_rank = comm.rank();
    int mpi_size = comm.size();
    int mpi_root = 0;

    atlas::vector<int> partition_global;
    ATLAS_TRACE_SCOPE("Gather global partition array") {
        Field field_partition_global = fs.createField(fs.partition(), option::global(mpi_root));
        fs.gather(fs.partition(), field_partition_global);
        int partition_global_size = field_partition_global.size();
        comm.broadcast(partition_global_size, mpi_root);

        partition_global.resize(partition_global_size);
        if (mpi_rank == mpi_root) {
            int* data = field_partition_global.array().host_data<int>();
            partition_global.assign(data, data+partition_global.size());
        }
        comm.broadcast(partition_global.data(), partition_global.size(), mpi_root);
    }

    atlas::grid::Distribution distribution(mpi_size, std::move(partition_global));
    locate(fs, distribution, size, global_index, partition, remote_index, remote_index_base);
}


inline void locate(const FunctionSpace& fs, const std::vector<gidx_t>& global_index, std::vector<int>& partition, std::vector<idx_t>& remote_index, idx_t remote_index_base) {
    auto size = global_index.size();
    partition.resize(size);
    remote_index.resize(size);
    locate(fs, size, global_index.data(), partition.data(), remote_index.data(), remote_index_base);
}


inline void locate(const FunctionSpace& fs, const grid::Distribution& distribution, const std::vector<gidx_t>& global_index, std::vector<int>& partition, std::vector<idx_t>& remote_index, idx_t remote_index_base) {
    auto size = global_index.size();
    partition.resize(size);
    remote_index.resize(size);
    locate(fs, distribution, size, global_index.data(), partition.data(), remote_index.data(), remote_index_base);
}

} // namespace atlas::util
