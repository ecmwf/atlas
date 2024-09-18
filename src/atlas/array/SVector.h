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

#include <algorithm>
#include <cassert>
#include <cstddef>

#include "atlas/library/config.h"
#include "atlas/util/Allocate.h"
#if ATLAS_HOST_COMPILE
#include "atlas/runtime/Exception.h"
#endif

namespace atlas {
namespace array {

//------------------------------------------------------------------------------

template <typename T>
class SVector {
public:
    ATLAS_HOST_DEVICE
    SVector(): data_(nullptr), size_(0), externally_allocated_(false) {}

    ATLAS_HOST_DEVICE
    SVector(const T* data, const idx_t size): data_(data), size_(size), externally_allocated_(true) {}

#if ATLAS_HIC_COMPILER
    // Note that this does not copy!!! It is mainly intended to be passed to a CUDA kernel which requires value semantics for this class
    ATLAS_HOST_DEVICE
    SVector(SVector const& other): data_(other.data_), size_(other.size_), externally_allocated_(true) {}
#endif

    ATLAS_HOST_DEVICE
    SVector(SVector&& other):
        data_(other.data_), size_(other.size_), externally_allocated_(other.externally_allocated_) {
        other.data_                 = nullptr;
        other.size_                 = 0;
        other.externally_allocated_ = true;
    }

    ATLAS_HOST_DEVICE
    SVector& operator=(SVector&& other) {
        data_                       = other.data_;
        size_                       = other.size_;
        externally_allocated_       = other.externally_allocated_;
        other.data_                 = nullptr;
        other.size_                 = 0;
        other.externally_allocated_ = true;
        return *this;
    }

    ATLAS_HOST_DEVICE
    SVector(T* data, idx_t size): data_(data), size_(size), externally_allocated_(true) {}

    SVector(idx_t N): data_(nullptr), size_(N), externally_allocated_(false) { allocate(data_, N); }

    ATLAS_HOST_DEVICE
    ~SVector() { clear(); }

    ATLAS_HOST_DEVICE
    void clear() {
        if (data_ && !externally_allocated_) {
#if ATLAS_HOST_COMPILE
            deallocate(data_, size_);
#endif
        }
        data_                 = nullptr;
        size_                 = 0;
        externally_allocated_ = false;
    }

    void insert(idx_t pos, idx_t dimsize) {
        T* data = nullptr;
        allocate(data, size_ + dimsize);
        for (idx_t c = 0; c < pos; ++c) {
            data[c] = std::move(data_[c]);
        }
        for (idx_t c = pos; c < size_; ++c) {
            data[c + dimsize] = std::move(data_[c]);
        }
        deallocate(data_, size_);
        data_ = data;
        size_ += dimsize;
    }

    size_t footprint() const { return sizeof(T) * size_; }

    ATLAS_HOST_DEVICE
    T* data() { return data_; }

    ATLAS_HOST_DEVICE
    T const* data() const { return data_; }

    ATLAS_HOST_DEVICE
    T& operator()(const idx_t idx) {
        //assert(data_ && idx < size_);
        return data_[idx];
    }
    ATLAS_HOST_DEVICE
    T const& operator()(const idx_t idx) const {
        //assert(data_ && idx < size_);
        return data_[idx];
    }

    ATLAS_HOST_DEVICE
    T& operator[](const idx_t idx) {
        //assert(data_ && idx < size_);
        return data_[idx];
    }
    ATLAS_HOST_DEVICE
    T const& operator[](const idx_t idx) const {
        //assert(data_ && idx < size_);
        return data_[idx];
    }

    ATLAS_HOST_DEVICE
    idx_t size() const { return size_; }

    void resize_impl(idx_t N) {
        if (N == size_)
            return;

        T* d_ = nullptr;
        allocate(d_, N);
        for (idx_t c = 0; c < std::min(size_, N); ++c) {
            d_[c] = std::move(data_[c]);
        }
        deallocate(data_, size_);
        data_ = d_;
    }

    void resize(idx_t N) {
#if ATLAS_HOST_COMPILE
        ATLAS_ASSERT(not externally_allocated_, "Cannot resize externally allocated (or wrapped) data");
#endif
        resize_impl(N);
        size_ = N;
    }

private:
    static void allocate(T*& ptr, idx_t size) {
        if (size > 0) {
            util::allocate_managedmem(ptr, size);
            for (idx_t c = 0; c < size; ++c) {
                new (ptr + c) T();
            }
        }
    }
    static void deallocate(T*& ptr, idx_t size) {
        if (ptr) {
            for (idx_t c = 0; c < size; ++c) {
                ptr[c].~T();
            }
            util::delete_managedmem(ptr, size);
        }
    }

private:
    T* data_;
    idx_t size_;
    bool externally_allocated_;
};

//------------------------------------------------------------------------------

}  // namespace array
}  // namespace atlas
