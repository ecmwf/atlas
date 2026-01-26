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

#include <stdexcept>
#include <string>
#include <vector>

#include "atlas/array/ArrayView.h"
#include "atlas/array/SVector.h"
#include "atlas/runtime/Exception.h"

namespace atlas::parallel::detail {

template <int ParallelDim, typename DATA_TYPE, int RANK>
struct HostPacker {
    static void pack(const int sendcnt, array::SVector<int> const& sendmap,
                     const array::ArrayView<DATA_TYPE, RANK>& array, DATA_TYPE* send_buffer, int send_buffer_size);
    static void unpack(const int recvcnt, array::SVector<int> const& recvmap, const DATA_TYPE* recv_buffer, int recv_buffer_size,
                       array::ArrayView<DATA_TYPE, RANK>& array);
};

template <int ParallelDim, typename DATA_TYPE, int RANK>
struct HostAdjointPacker {
    static void unpack(const int recvcnt, array::SVector<int> const& recvmap, const DATA_TYPE* recv_buffer,
                       int recv_buffer_size, array::ArrayView<DATA_TYPE, RANK>& array);
};

template <int ParallelDim, typename DATA_TYPE, int RANK>
struct HostZeroer {
    static void zero(const int sendcnt, array::SVector<int> const& sendmap, array::ArrayView<DATA_TYPE, RANK>& array,
                       DATA_TYPE* recv_buffer, int recv_buffer_size);
};

#if ATLAS_HAVE_GPU
template <int ParallelDim, typename DATA_TYPE, int RANK>
struct DevicePacker {
    static void pack(const int sendcnt, array::SVector<int> const& sendmap,
                     const array::ArrayView<DATA_TYPE, RANK>& array,
                     DATA_TYPE* send_buffer, int send_buffer_size);

    static void unpack(const int sendcnt, array::SVector<int> const& recvmap,
                       const DATA_TYPE* recv_buffer, int recv_buffer_size,
                       array::ArrayView<DATA_TYPE, RANK>& array);
};
#else
template <int ParallelDim, typename DATA_TYPE, int RANK>
using DevicePacker = HostPacker<ParallelDim, DATA_TYPE, RANK>;
#endif



template <int ParallelDim, typename DATA_TYPE, int RANK>
struct Packer {
    static void pack(const int sendcnt, array::SVector<int> const& sendmap,
                     const array::ArrayView<DATA_TYPE, RANK>& array,
                     DATA_TYPE* send_buffer, int send_buffer_size,
                     bool on_device) {
        if (not on_device) {
            HostPacker<ParallelDim, DATA_TYPE, RANK>::pack(sendcnt, sendmap, array, send_buffer, send_buffer_size);
        }
        else {
            DevicePacker<ParallelDim, DATA_TYPE, RANK>::pack(sendcnt, sendmap, array, send_buffer, send_buffer_size);
        }
    }

    static void unpack(const int recvcnt, array::SVector<int> const& recvmap,
                       const DATA_TYPE* recv_buffer, int recv_buffer_size,
                       array::ArrayView<DATA_TYPE, RANK>& array,
                       bool on_device) {
        if (not on_device) {
            HostPacker<ParallelDim, DATA_TYPE, RANK>::unpack(recvcnt, recvmap, recv_buffer, recv_buffer_size, array);
        }
        else {
            DevicePacker<ParallelDim, DATA_TYPE, RANK>::unpack(recvcnt, recvmap, recv_buffer, recv_buffer_size, array);
        }
    }
};

template <int ParallelDim, typename DATA_TYPE, int RANK>
struct AdjointPacker {
    static void pack(const int sendcnt, array::SVector<int> const& sendmap,
                     const array::ArrayView<DATA_TYPE, RANK>& array,
                     DATA_TYPE* send_buffer, int send_buffer_size,
                     bool on_device) {
        Packer<ParallelDim, DATA_TYPE, RANK>::pack(sendcnt, sendmap, array, send_buffer, send_buffer_size, on_device);
    }
    static void unpack(const int recvcnt, array::SVector<int> const& recvmap, const DATA_TYPE* recv_buffer,
                       int recv_buffer_size, array::ArrayView<DATA_TYPE, RANK>& array,
                       bool on_device) {
        if (on_device) {
            ATLAS_NOTIMPLEMENTED;
        }
        else {
            HostAdjointPacker<ParallelDim, DATA_TYPE, RANK>::unpack(recvcnt, recvmap, recv_buffer, recv_buffer_size, array);
        }
    }
};


template <int ParallelDim, typename DATA_TYPE, int RANK>
struct Zeroer {
    static void zero(const int sendcnt, array::SVector<int> const& sendmap, array::ArrayView<DATA_TYPE, RANK>& array,
                       DATA_TYPE* recv_buffer, int recv_buffer_size,
                       bool on_device) {
        if (on_device) {
            ATLAS_NOTIMPLEMENTED;
        }
        else {
            HostZeroer<ParallelDim, DATA_TYPE, RANK>::zero(sendcnt, sendmap, array, recv_buffer, recv_buffer_size);
        }
    }
};

//----------------------------------------------------------------------------------------------------------------------

}  // namespace atlas::parallel::detail
