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

#include <type_traits>
#include <utility>  // std::swap

#include "atlas/library/config.h"
#include "atlas/parallel/omp/copy.h"
#include "atlas/parallel/omp/fill.h"
#include "atlas/runtime/Exception.h"

namespace atlas {

template <typename T>
class vector {
public:
    using value_type     = T;
    using iterator       = T*;
    using const_iterator = T const*;

public:
    vector() = default;

    template <typename size_t, typename std::enable_if<std::is_integral<size_t>::value, int>::type = 0>
    vector(size_t size) {
        resize(size);
    }

    template <typename size_t, typename std::enable_if<std::is_integral<size_t>::value, int>::type = 0>
    vector(size_t size, const value_type& value): vector(size) {
        assign(size, value);
    }

    vector(const vector& other) { assign(other.data_, other.data_ + other.size_); }

    vector(vector&& other) {
        std::swap(data_, other.data_);
        std::swap(size_, other.size_);
        std::swap(capacity_, other.capacity_);
    }

    vector& operator=(vector other) {
        std::swap(data_, other.data_);
        std::swap(size_, other.size_);
        std::swap(capacity_, other.capacity_);
        return *this;
    }

    template <typename T2>
    vector(const std::initializer_list<T2>& list) {
        assign(list.begin(), list.end());
    }

    ~vector() {
        if (data_) {
            delete[] data_;
        }
    }

    template <typename idx_t, typename std::enable_if<std::is_integral<idx_t>::value, int>::type = 0>
    T& at(idx_t i) noexcept(false) {
        if (i >= size_) {
            throw_OutOfRange("atlas::vector", i, size_);
        }
        return data_[i];
    }

    template <typename idx_t, typename std::enable_if<std::is_integral<idx_t>::value, int>::type = 0>
    T const& at(idx_t i) const noexcept(false) {
        if (i >= size_) {
            throw_OutOfRange("atlas::vector", i, size_);
        }
        return data_[i];
    }

    template <typename idx_t, typename std::enable_if<std::is_integral<idx_t>::value, int>::type = 0>
    T& operator[](idx_t i) {
#if ATLAS_VECTOR_BOUNDS_CHECKING
        return at(i);
#else
        return data_[i];
#endif
    }

    template <typename idx_t, typename std::enable_if<std::is_integral<idx_t>::value, int>::type = 0>
    T const& operator[](idx_t i) const {
#if ATLAS_VECTOR_BOUNDS_CHECKING
        return at(i);
#else
        return data_[i];
#endif
    }

    const T* data() const { return data_; }

    T* data() { return data_; }

    idx_t size() const { return size_; }

    idx_t capacity() const { return capacity_; }

    template <typename Size, typename std::enable_if<std::is_integral<Size>::value, int>::type = 0>
    void assign(Size n, const value_type& value) {
        resize(n);
        omp::fill(begin(), begin() + n, value);
    }

    template <typename Iter, typename std::enable_if<!std::is_integral<Iter>::value, int>::type = 0>
    void assign(const Iter& first, const Iter& last) {
        size_t size = std::distance(first, last);
        resize(size);
        omp::copy(first, last, begin());
    }

    template <typename Size, typename std::enable_if<std::is_integral<Size>::value, int>::type = 0>
    void reserve(Size size) {
        if (capacity_ != 0)
            ATLAS_NOTIMPLEMENTED;
        data_     = new T[size];
        capacity_ = size;
    }
    template <typename size_t, typename std::enable_if<std::is_integral<size_t>::value, int>::type = 0>
    void resize(size_t size) {
        if (static_cast<idx_t>(size) > 0) {
            if (capacity_ == 0) {
                reserve(size);
            }
            if (static_cast<idx_t>(size) > capacity_) {
                ATLAS_NOTIMPLEMENTED;
            }
            size_ = size;
        }
    }
    const_iterator begin() const { return data_; }
    const_iterator end() const { return data_ + size_; }
    iterator begin() { return data_; }
    iterator end() { return data_ + size_; }
    const_iterator cbegin() const { return data_; }
    const_iterator cend() const { return data_ + size_; }

private:
    value_type* data_{nullptr};
    idx_t size_{0};
    idx_t capacity_{0};
};

}  // namespace atlas
