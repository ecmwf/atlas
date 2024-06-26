/*
* (C) Copyright 2013 ECMWF.
*
* This software is licensed under the terms of the Apache Licence Version 2.0
* which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
* In applying this licence, ECMWF does not waive the privileges and immunities
* granted to it by virtue of its status as an intergovernmental organisation nor
* does it submit to any jurisdiction.
*/

#pragma once

#include <cassert>
#include <cstddef>

#include "atlas/library/config.h"
#include "atlas/runtime/Exception.h"

#if ATLAS_HAVE_CUDA
#include <cuda_runtime.h>
#endif

namespace atlas {
namespace array {

//------------------------------------------------------------------------------

template <typename T>
class Vector {
public:
    Vector(idx_t N = 0): data_(N ? new T[N] : nullptr), data_gpu_(nullptr), size_(N) {}

    void resize_impl(idx_t N) {
        ATLAS_ASSERT(not data_gpu_, "we can not resize a vector after has been cloned to device");
        ATLAS_ASSERT(N >= size_);
        if (N == size_)
            return;

        T* d_ = new T[N];
        for (idx_t c = 0; c < size_; ++c) {
            d_[c] = data_[c];
        }
        if (data_)
            delete[] data_;
        data_ = d_;
    }

    void resize(idx_t N) {
        resize_impl(N);
        size_ = N;
    }

    void resize(idx_t N, T val) {
        resize_impl(N);
        for (idx_t c = size_; c < N; ++c) {
            data_[c] = val;
        }

        size_ = N;
    }

    void updateDevice() {
        if (!data_gpu_) {
#if ATLAS_HAVE_CUDA
            ::cudaMalloc((void**)(&data_gpu_), sizeof(T*) * size_);

            T* buff = new T[size_];

            for (idx_t i = 0; i < size(); ++i) {
                data_[i]->updateDevice();
                buff[i] = data_[i]->gpu_object_ptr();
            }
            ::cudaMemcpy(data_gpu_, buff, sizeof(T*) * size_, cudaMemcpyHostToDevice);
            delete buff;
#else
            data_gpu_ = data_;
#endif
            size_gpu_ = size_;
        }
        else {
            ATLAS_ASSERT(size_gpu_ == size_);
#if ATLAS_HAVE_CUDA
            for (idx_t i = 0; i < size(); ++i) {
                data_[i]->updateDevice();
                assert(data_gpu_[i] == data_[i]->gpu_object_ptr());
            }
#endif
        }
    }
    void updateHost() {
        ATLAS_ASSERT(data_gpu_ != nullptr);

#if ATLAS_HAVE_CUDA

        for (idx_t i = 0; i < size(); ++i) {
            data_[i]->updateHost();
        }
#endif
    }

    T* gpu_object_ptr() { return data_gpu_; }

    T* data() { return data_; }
    T* data_gpu() { return data_gpu_; }

    idx_t size() const { return size_; }

private:
    T* data_;
    T* data_gpu_;
    idx_t size_;
    idx_t size_gpu_;
};

template <typename T>
class VectorView {
public:
    VectorView(): vector_(nullptr), data_(nullptr), size_(0) {}
    VectorView(Vector<T> const& vector, T* data): vector_(&vector), data_(data), size_(vector.size()) {}

    ATLAS_HOST_DEVICE
    T& operator[](idx_t idx) {
        assert(idx < size_);

        return data_[idx];
    }

    ATLAS_HOST_DEVICE
    T const& operator[](idx_t idx) const {
        assert(idx < size_);
        return data_[idx];
    }
    ATLAS_HOST_DEVICE
    T base() { return *data_; }

    ATLAS_HOST_DEVICE
    idx_t size() const { return size_; }

    bool is_valid(Vector<T>& vector) { return (&vector) == vector_ && (data_ != nullptr); }

public:
    Vector<T> const* vector_;
    T* data_;
    idx_t size_;
};

template <typename T>
VectorView<T> make_host_vector_view(Vector<T> vector_) {
    return VectorView<T>(vector_, vector_.data());
}

template <typename T>
VectorView<T> make_device_vector_view(Vector<T> vector_) {
    return VectorView<T>(vector_, vector_.data_gpu());
}

//------------------------------------------------------------------------------

}  // namespace array
}  // namespace atlas
