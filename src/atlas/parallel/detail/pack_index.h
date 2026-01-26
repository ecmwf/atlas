/*
 * (C) Copyright 2025- ECMWF.
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

namespace atlas::parallel::detail {

template <int ParallelDim, int Cnt, int CurrentDim>
struct pack_index {
    template <typename DATA_TYPE, int RANK, typename... Idx>
    ATLAS_HOST_DEVICE static void apply(idx_t& buf_idx, const idx_t node_idx,
                                        const array::ArrayView<DATA_TYPE, RANK>& field, DATA_TYPE* send_buffer,
                                        Idx... idxs) {
        for (idx_t i = 0; i < field.template shape<CurrentDim>(); ++i) {
            pack_index<ParallelDim, Cnt - 1, CurrentDim + 1>::apply(buf_idx, node_idx, field, send_buffer,
                                                                          idxs..., i);
        }
    }
};

template <int ParallelDim>
struct pack_index<ParallelDim, 0, ParallelDim> {
    template <typename DATA_TYPE, int RANK, typename... Idx>
    ATLAS_HOST_DEVICE static void apply(idx_t& buf_idx, const idx_t node_idx,
                                        const array::ArrayView<DATA_TYPE, RANK>& field, DATA_TYPE* send_buffer,
                                        Idx... idxs) {
        send_buffer[buf_idx++] = field(idxs...);
    }
};

template <int ParallelDim, int Cnt>
struct pack_index<ParallelDim, Cnt, ParallelDim> {
    template <typename DATA_TYPE, int RANK, typename... Idx>
    ATLAS_HOST_DEVICE static void apply(idx_t& buf_idx, const idx_t node_idx,
                                        const array::ArrayView<DATA_TYPE, RANK>& field, DATA_TYPE* send_buffer,
                                        Idx... idxs) {
        pack_index<ParallelDim, Cnt - 1, ParallelDim + 1>::apply(buf_idx, node_idx, field, send_buffer, idxs...,
                                                                       node_idx);
    }
};

template <int ParallelDim, int CurrentDim>
struct pack_index<ParallelDim, 0, CurrentDim> {
    template <typename DATA_TYPE, int RANK, typename... Idx>
    ATLAS_HOST_DEVICE static void apply(idx_t& buf_idx, const idx_t node_idx,
                                        const array::ArrayView<DATA_TYPE, RANK>& field, DATA_TYPE* send_buffer,
                                        Idx... idxs) {
        send_buffer[buf_idx++] = field(idxs...);
    }
};

template <int ParallelDim, int Cnt, int CurrentDim>
struct unpack_index {
    template <typename DATA_TYPE, int RANK, typename... Idx>
    ATLAS_HOST_DEVICE static void apply(idx_t& buf_idx, const idx_t node_idx, const DATA_TYPE* recv_buffer,
                                        array::ArrayView<DATA_TYPE, RANK>& field, Idx... idxs) {
        for (idx_t i = 0; i < field.template shape<CurrentDim>(); ++i) {
            unpack_index<ParallelDim, Cnt - 1, CurrentDim + 1>::apply(buf_idx, node_idx, recv_buffer, field,
                                                                            idxs..., i);
        }
    }
};

template <int ParallelDim>
struct unpack_index<ParallelDim, 0, ParallelDim> {
    template <typename DATA_TYPE, int RANK, typename... Idx>
    ATLAS_HOST_DEVICE static void apply(idx_t& buf_idx, const idx_t node_idx, const DATA_TYPE* recv_buffer,
                                        array::ArrayView<DATA_TYPE, RANK>& field, Idx... idxs) {
        field(idxs...) = recv_buffer[buf_idx++];
    }
};

template <int ParallelDim, int Cnt>
struct unpack_index<ParallelDim, Cnt, ParallelDim> {
    template <typename DATA_TYPE, int RANK, typename... Idx>
    ATLAS_HOST_DEVICE static void apply(idx_t& buf_idx, const idx_t node_idx, const DATA_TYPE* recv_buffer,
                                        array::ArrayView<DATA_TYPE, RANK>& field, Idx... idxs) {
        unpack_index<ParallelDim, Cnt - 1, ParallelDim + 1>::apply(buf_idx, node_idx, recv_buffer, field, idxs...,
                                                                         node_idx);
    }
};

template <int ParallelDim, int CurrentDim>
struct unpack_index<ParallelDim, 0, CurrentDim> {
    template <typename DATA_TYPE, int RANK, typename... Idx>
    ATLAS_HOST_DEVICE static void apply(idx_t& buf_idx, const idx_t node_idx, const DATA_TYPE* recv_buffer,
                                        array::ArrayView<DATA_TYPE, RANK>& field, Idx... idxs) {
        field(idxs...) = recv_buffer[buf_idx++];
    }
};

}  // namespace atlas::parallel::detail
