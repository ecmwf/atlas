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
struct adjoint_unpack_index {
    template <typename DATA_TYPE, int RANK, typename... Idx>
    ATLAS_HOST_DEVICE static void apply(idx_t& buf_idx, const idx_t node_idx, const DATA_TYPE* recv_buffer,
                                        array::ArrayView<DATA_TYPE, RANK>& field, Idx... idxs) {
        for (idx_t i = 0; i < field.template shape<CurrentDim>(); ++i) {
            adjoint_unpack_index<ParallelDim, Cnt - 1, CurrentDim + 1>::apply(buf_idx, node_idx, recv_buffer,
                                                                                    field, idxs..., i);
        }
    }
};


template <int ParallelDim>
struct adjoint_unpack_index<ParallelDim, 0, ParallelDim> {
    template <typename DATA_TYPE, int RANK, typename... Idx>
    ATLAS_HOST_DEVICE static void apply(idx_t& buf_idx, const idx_t node_idx, const DATA_TYPE* recv_buffer,
                                        array::ArrayView<DATA_TYPE, RANK>& field, Idx... idxs) {
        field(idxs...) += recv_buffer[buf_idx++];
    }
};

template <int ParallelDim, int Cnt>
struct adjoint_unpack_index<ParallelDim, Cnt, ParallelDim> {
    template <typename DATA_TYPE, int RANK, typename... Idx>
    ATLAS_HOST_DEVICE static void apply(idx_t& buf_idx, const idx_t node_idx, const DATA_TYPE* recv_buffer,
                                        array::ArrayView<DATA_TYPE, RANK>& field, Idx... idxs) {
        adjoint_unpack_index<ParallelDim, Cnt - 1, ParallelDim + 1>::apply(buf_idx, node_idx, recv_buffer, field,
                                                                                 idxs..., node_idx);
    }
};

template <int ParallelDim, int CurrentDim>
struct adjoint_unpack_index<ParallelDim, 0, CurrentDim> {
    template <typename DATA_TYPE, int RANK, typename... Idx>
    ATLAS_HOST_DEVICE static void apply(idx_t& buf_idx, const idx_t node_idx, const DATA_TYPE* recv_buffer,
                                        array::ArrayView<DATA_TYPE, RANK>& field, Idx... idxs) {
        field(idxs...) += recv_buffer[buf_idx++];
    }
};

}  // namespace atlas::parallel::detail
