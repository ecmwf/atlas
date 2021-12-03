/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#pragma once

#include "atlas/array/ArrayView.h"
#include "atlas/array/SVector.h"

namespace atlas {
namespace parallel {

template <int ParallelDim, int Cnt, int CurrentDim>
struct halo_packer_impl {
    template <typename DATA_TYPE, int RANK, typename... Idx>
    ATLAS_HOST_DEVICE static void apply(idx_t& buf_idx, const idx_t node_idx,
                                        const array::ArrayView<DATA_TYPE, RANK>& field, DATA_TYPE* send_buffer,
                                        Idx... idxs) {
        for (idx_t i = 0; i < field.template shape<CurrentDim>(); ++i) {
            halo_packer_impl<ParallelDim, Cnt - 1, CurrentDim + 1>::apply(buf_idx, node_idx, field, send_buffer,
                                                                          idxs..., i);
        }
    }
};

template <int ParallelDim>
struct halo_packer_impl<ParallelDim, 0, ParallelDim> {
    template <typename DATA_TYPE, int RANK, typename... Idx>
    ATLAS_HOST_DEVICE static void apply(idx_t& buf_idx, const idx_t node_idx,
                                        const array::ArrayView<DATA_TYPE, RANK>& field, DATA_TYPE* send_buffer,
                                        Idx... idxs) {
        send_buffer[buf_idx++] = field(idxs...);
    }
};

template <int ParallelDim, int Cnt>
struct halo_packer_impl<ParallelDim, Cnt, ParallelDim> {
    template <typename DATA_TYPE, int RANK, typename... Idx>
    ATLAS_HOST_DEVICE static void apply(idx_t& buf_idx, const idx_t node_idx,
                                        const array::ArrayView<DATA_TYPE, RANK>& field, DATA_TYPE* send_buffer,
                                        Idx... idxs) {
        halo_packer_impl<ParallelDim, Cnt - 1, ParallelDim + 1>::apply(buf_idx, node_idx, field, send_buffer, idxs...,
                                                                       node_idx);
    }
};

template <int ParallelDim, int CurrentDim>
struct halo_packer_impl<ParallelDim, 0, CurrentDim> {
    template <typename DATA_TYPE, int RANK, typename... Idx>
    ATLAS_HOST_DEVICE static void apply(idx_t& buf_idx, const idx_t node_idx,
                                        const array::ArrayView<DATA_TYPE, RANK>& field, DATA_TYPE* send_buffer,
                                        Idx... idxs) {
        send_buffer[buf_idx++] = field(idxs...);
    }
};

template <int ParallelDim, int Cnt, int CurrentDim>
struct halo_unpacker_impl {
    template <typename DATA_TYPE, int RANK, typename... Idx>
    ATLAS_HOST_DEVICE static void apply(idx_t& buf_idx, const idx_t node_idx, const DATA_TYPE* recv_buffer,
                                        array::ArrayView<DATA_TYPE, RANK>& field, Idx... idxs) {
        for (idx_t i = 0; i < field.template shape<CurrentDim>(); ++i) {
            halo_unpacker_impl<ParallelDim, Cnt - 1, CurrentDim + 1>::apply(buf_idx, node_idx, recv_buffer, field,
                                                                            idxs..., i);
        }
    }
};

template <int ParallelDim>
struct halo_unpacker_impl<ParallelDim, 0, ParallelDim> {
    template <typename DATA_TYPE, int RANK, typename... Idx>
    ATLAS_HOST_DEVICE static void apply(idx_t& buf_idx, const idx_t node_idx, const DATA_TYPE* recv_buffer,
                                        array::ArrayView<DATA_TYPE, RANK>& field, Idx... idxs) {
        field(idxs...) = recv_buffer[buf_idx++];
    }
};

template <int ParallelDim, int Cnt>
struct halo_unpacker_impl<ParallelDim, Cnt, ParallelDim> {
    template <typename DATA_TYPE, int RANK, typename... Idx>
    ATLAS_HOST_DEVICE static void apply(idx_t& buf_idx, const idx_t node_idx, const DATA_TYPE* recv_buffer,
                                        array::ArrayView<DATA_TYPE, RANK>& field, Idx... idxs) {
        halo_unpacker_impl<ParallelDim, Cnt - 1, ParallelDim + 1>::apply(buf_idx, node_idx, recv_buffer, field, idxs...,
                                                                         node_idx);
    }
};

template <int ParallelDim, int CurrentDim>
struct halo_unpacker_impl<ParallelDim, 0, CurrentDim> {
    template <typename DATA_TYPE, int RANK, typename... Idx>
    ATLAS_HOST_DEVICE static void apply(idx_t& buf_idx, const idx_t node_idx, const DATA_TYPE* recv_buffer,
                                        array::ArrayView<DATA_TYPE, RANK>& field, Idx... idxs) {
        field(idxs...) = recv_buffer[buf_idx++];
    }
};

}  // namespace parallel
}  // namespace atlas
