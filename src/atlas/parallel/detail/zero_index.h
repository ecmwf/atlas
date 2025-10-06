/*
 * (C) British Crown Copyright 2020 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include "atlas/array/ArrayView.h"
#include "atlas/array/SVector.h"

namespace atlas::parallel::detail {

template <int ParallelDim, int Cnt, int CurrentDim>
struct zero_index {
    template <typename DATA_TYPE, int RANK, typename... Idx>
    ATLAS_HOST_DEVICE static void apply(idx_t& buf_idx, const idx_t node_idx, array::ArrayView<DATA_TYPE, RANK>& field,
                                        DATA_TYPE* recv_buffer, Idx... idxs) {
        for (idx_t i = 0; i < field.template shape<CurrentDim>(); ++i) {
            zero_index<ParallelDim, Cnt - 1, CurrentDim + 1>::apply(buf_idx, node_idx, field, recv_buffer,
                                                                          idxs..., i);
        }
    }
};

template <int ParallelDim>
struct zero_index<ParallelDim, 0, ParallelDim> {
    template <typename DATA_TYPE, int RANK, typename... Idx>
    ATLAS_HOST_DEVICE static void apply(idx_t& buf_idx, const idx_t node_idx, array::ArrayView<DATA_TYPE, RANK>& field,
                                        DATA_TYPE* recv_buffer, Idx... idxs) {
        field(idxs...) = 0;
    }
};

template <int ParallelDim, int Cnt>
struct zero_index<ParallelDim, Cnt, ParallelDim> {
    template <typename DATA_TYPE, int RANK, typename... Idx>
    ATLAS_HOST_DEVICE static void apply(idx_t& buf_idx, const idx_t node_idx, array::ArrayView<DATA_TYPE, RANK>& field,
                                        DATA_TYPE* recv_buffer, Idx... idxs) {
        zero_index<ParallelDim, Cnt - 1, ParallelDim + 1>::apply(buf_idx, node_idx, field, recv_buffer, idxs...,
                                                                       node_idx);
    }
};

template <int ParallelDim, int CurrentDim>
struct zero_index<ParallelDim, 0, CurrentDim> {
    template <typename DATA_TYPE, int RANK, typename... Idx>
    ATLAS_HOST_DEVICE static void apply(idx_t& buf_idx, const idx_t node_idx, array::ArrayView<DATA_TYPE, RANK>& field,
                                        DATA_TYPE* recv_buffer, Idx... idxs) {
        field(idxs...) = 0;
    }
};

}  // namespace atlas::parallel::detail
