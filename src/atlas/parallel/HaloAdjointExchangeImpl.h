/*
 * (C) British Crown Copyright 2020 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include "atlas/array/ArrayView.h"
#include "atlas/array/SVector.h"

namespace atlas {
namespace parallel {

template <int ParallelDim, int Cnt, int CurrentDim>
struct halo_zeroer_impl {
    template <typename DATA_TYPE, int RANK, typename... Idx>
    ATLAS_HOST_DEVICE static void apply(idx_t& buf_idx, const idx_t node_idx, array::ArrayView<DATA_TYPE, RANK>& field,
                                        DATA_TYPE* recv_buffer, Idx... idxs) {
        for (idx_t i = 0; i < field.template shape<CurrentDim>(); ++i) {
            halo_zeroer_impl<ParallelDim, Cnt - 1, CurrentDim + 1>::apply(buf_idx, node_idx, field, recv_buffer,
                                                                          idxs..., i);
        }
    }
};

template <int ParallelDim>
struct halo_zeroer_impl<ParallelDim, 0, ParallelDim> {
    template <typename DATA_TYPE, int RANK, typename... Idx>
    ATLAS_HOST_DEVICE static void apply(idx_t& buf_idx, const idx_t node_idx, array::ArrayView<DATA_TYPE, RANK>& field,
                                        DATA_TYPE* recv_buffer, Idx... idxs) {
        field(idxs...) = 0;
    }
};

template <int ParallelDim, int Cnt>
struct halo_zeroer_impl<ParallelDim, Cnt, ParallelDim> {
    template <typename DATA_TYPE, int RANK, typename... Idx>
    ATLAS_HOST_DEVICE static void apply(idx_t& buf_idx, const idx_t node_idx, array::ArrayView<DATA_TYPE, RANK>& field,
                                        DATA_TYPE* recv_buffer, Idx... idxs) {
        halo_zeroer_impl<ParallelDim, Cnt - 1, ParallelDim + 1>::apply(buf_idx, node_idx, field, recv_buffer, idxs...,
                                                                       node_idx);
    }
};

template <int ParallelDim, int CurrentDim>
struct halo_zeroer_impl<ParallelDim, 0, CurrentDim> {
    template <typename DATA_TYPE, int RANK, typename... Idx>
    ATLAS_HOST_DEVICE static void apply(idx_t& buf_idx, const idx_t node_idx, array::ArrayView<DATA_TYPE, RANK>& field,
                                        DATA_TYPE* recv_buffer, Idx... idxs) {
        field(idxs...) = 0;
    }
};

template <int ParallelDim, int Cnt, int CurrentDim>
struct halo_adjoint_unpacker_impl {
    template <typename DATA_TYPE, int RANK, typename... Idx>
    ATLAS_HOST_DEVICE static void apply(idx_t& buf_idx, const idx_t node_idx, const DATA_TYPE* recv_buffer,
                                        array::ArrayView<DATA_TYPE, RANK>& field, Idx... idxs) {
        for (idx_t i = 0; i < field.template shape<CurrentDim>(); ++i) {
            halo_adjoint_unpacker_impl<ParallelDim, Cnt - 1, CurrentDim + 1>::apply(buf_idx, node_idx, recv_buffer,
                                                                                    field, idxs..., i);
        }
    }
};


template <int ParallelDim>
struct halo_adjoint_unpacker_impl<ParallelDim, 0, ParallelDim> {
    template <typename DATA_TYPE, int RANK, typename... Idx>
    ATLAS_HOST_DEVICE static void apply(idx_t& buf_idx, const idx_t node_idx, const DATA_TYPE* recv_buffer,
                                        array::ArrayView<DATA_TYPE, RANK>& field, Idx... idxs) {
        field(idxs...) += recv_buffer[buf_idx++];
    }
};

template <int ParallelDim, int Cnt>
struct halo_adjoint_unpacker_impl<ParallelDim, Cnt, ParallelDim> {
    template <typename DATA_TYPE, int RANK, typename... Idx>
    ATLAS_HOST_DEVICE static void apply(idx_t& buf_idx, const idx_t node_idx, const DATA_TYPE* recv_buffer,
                                        array::ArrayView<DATA_TYPE, RANK>& field, Idx... idxs) {
        halo_adjoint_unpacker_impl<ParallelDim, Cnt - 1, ParallelDim + 1>::apply(buf_idx, node_idx, recv_buffer, field,
                                                                                 idxs..., node_idx);
    }
};

template <int ParallelDim, int CurrentDim>
struct halo_adjoint_unpacker_impl<ParallelDim, 0, CurrentDim> {
    template <typename DATA_TYPE, int RANK, typename... Idx>
    ATLAS_HOST_DEVICE static void apply(idx_t& buf_idx, const idx_t node_idx, const DATA_TYPE* recv_buffer,
                                        array::ArrayView<DATA_TYPE, RANK>& field, Idx... idxs) {
        field(idxs...) += recv_buffer[buf_idx++];
    }
};

}  // namespace parallel
}  // namespace atlas
