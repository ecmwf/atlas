/*
 * (C) Copyright 2025 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/functionspace/Locate.h"

#include "atlas/functionspace/FunctionSpace.h"
#include "atlas/array/ArrayView.h"
#include "atlas/array/Array.h"
#include "atlas/option.h"
#include "atlas/parallel/mpi/mpi.h"

namespace atlas::functionspace {

// Figure out the global distribution using the functionspace partition fields from each rank.
static atlas::vector<int> allgather_distribution_array(const FunctionSpace& fs) {
    auto& comm = atlas::mpi::comm(fs.mpi_comm());
    int mpi_rank = comm.rank();
    int mpi_root = 0;

    atlas::vector<int> distribution_array;
    ATLAS_TRACE_SCOPE("Gather global partition array") {
        Field field_partition_global = fs.createField(fs.partition(), option::global(mpi_root));
        fs.gather(fs.partition(), field_partition_global);
        std::size_t partition_global_size = field_partition_global.size();
        comm.broadcast(partition_global_size, mpi_root);
        distribution_array.resize(partition_global_size);
        if (mpi_rank == mpi_root) {
            int* data = field_partition_global.array().host_data<int>();
            distribution_array.assign(data, data+partition_global_size);
        }
        comm.broadcast(distribution_array.data(), distribution_array.size(), mpi_root);
    }
    return distribution_array;
}

Locator::Locator(const FunctionSpace& fs) :
    fs_global_index_{fs.global_index()},
    fs_ghost_{fs.ghost()} {
    distribution_array_ = allgather_distribution_array(fs);
    distribution_function_ = [this](gidx_t j){ return distribution_array_[j]; };
    distribution_size_ = distribution_array_.size();
    mpi_comm_ = fs.mpi_comm();
}

void Locator::locate(
    // input
    span<const gidx_t> global_index, const gidx_t global_index_base,
    // output
    span<int> partition, const int partition_base,
    span<idx_t> remote_index, const idx_t remote_index_base) const {
    constexpr gidx_t fs_global_index_base_ = 1;
    constexpr int distribution_base_ = 0;
    if (distribution_array_.size()) {
        ::atlas::parallel::Locator::locate(
            // context
            mpi_comm_,
            array::make_view<const gidx_t,1>(fs_global_index_).as_mdspan(), fs_global_index_base_,
            array::make_view<const int,1>(fs_ghost_).as_mdspan(),
            span<const int>{distribution_array_.data(), distribution_array_.size()}, distribution_base_,
            // input
            global_index, global_index_base,
            // output
            partition, partition_base,
            remote_index, remote_index_base);
    }
    else {
        ::atlas::parallel::Locator::locate(
            // context
            mpi_comm_,
            array::make_view<const gidx_t,1>(fs_global_index_).as_mdspan(),
            fs_global_index_base_,
            array::make_view<const int,1>(fs_ghost_).as_mdspan(),
            fspan<const int>{&distribution_function_, distribution_size_}, distribution_base_,
            // input
            global_index, global_index_base,
            // output
            partition, partition_base,
            remote_index, remote_index_base);
    }
}

} // namespace atlas::functionspace
