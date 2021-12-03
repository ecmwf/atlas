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

#include "atlas/array/LocalView.h"
#include "atlas/parallel/mpi/mpi.h"

namespace atlas {
namespace mpi {

/// @brief Buffer<DATA_TYPE,SHAPE>
///
/// Class that keeps allocation of a MPI buffer including
/// counts and displacements, but with added index operator[]
/// that returns an array::ArrayView<DATA_TYPE,SHAPE> of the part
/// of the buffer for a processor index.
template <typename DATA_TYPE, int SHAPE>
struct Buffer;

// ----------------------------------------------------------------------------------

template <typename DATA_TYPE, int SHAPE>
struct Buffer : eckit::mpi::Buffer<DATA_TYPE> {};

template <typename DATA_TYPE>
class BufferView {
public:
    BufferView(DATA_TYPE* data, size_t size): data_(data), size_(size) {}
    size_t size() const { return size_; }
    const DATA_TYPE& operator()(const size_t i) const { return data_[i]; }
    const DATA_TYPE& operator[](const size_t i) const { return data_[i]; }

private:
    DATA_TYPE* data_;
    size_t size_;
};

template <typename DATA_TYPE>
struct Buffer<DATA_TYPE, 1> : public eckit::mpi::Buffer<DATA_TYPE> {
    Buffer(size_t size): eckit::mpi::Buffer<DATA_TYPE>(size) {}

    //                                   array::make_shape(
    //                                   eckit::mpi::Buffer<DATA_TYPE>::counts[p]
    //                                   ) );
    // }
    BufferView<DATA_TYPE> operator[](int p) {
        return BufferView<DATA_TYPE>(
            eckit::mpi::Buffer<DATA_TYPE>::buffer.data() + eckit::mpi::Buffer<DATA_TYPE>::displs[p],
            eckit::mpi::Buffer<DATA_TYPE>::counts[p]);
    }
};

// ----------------------------------------------------------------------------------

}  // namespace mpi
}  // namespace atlas
