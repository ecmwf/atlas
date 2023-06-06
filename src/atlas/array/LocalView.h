/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

/// @file LocalView.h
/// This file contains the LocalView class, a class that allows to wrap any
/// contiguous raw data into
/// a view which is accessible with multiple indices.

#pragma once

#include <cstddef>
#include <type_traits>
#include <array>

#include "atlas/array/ArrayDataStore.h"
#include "atlas/array/ArrayViewDefs.h"
#include "atlas/array/helpers/ArraySlicer.h"
#include "atlas/library/config.h"

//------------------------------------------------------------------------------------------------------

namespace atlas {
namespace array {

/// @brief Multi-dimensional access existing POD array pointer.
///
/// A LocalView is a wrapper around data that enables multidimensional access, and has the exact
/// same API as ArrayView.
///
/// The data may be strided.
///
/// ### Example 1:
///
/// @code{.cpp}
///     int[] array = { 1, 2, 3, 4, 5, 6, 7, 8, 9};
///     int[2] strides = { 3, 1 };
///     int[2] shape = { 3, 3 };
///     LocalView<int,2> matrix( array, shape, strides );
///     for( idx_t i=0; i<matrix.shape(0); ++i ) {
///         for( idx_t j=0; j<matrix.shape(1); ++j ) {
///             matrix(i,j) *= 10;
///         }
///     }
/// @endcode
///
/// Strides can also be omitted as for most common cases it can be inferred
/// from the shape.
///
/// ### Example 2:
///
/// @code{.cpp}
///     int[] array = { 1, 2, 3, 4, 5, 6, 7, 8, 9};
///     int[2] shape = { 3, 3 };
///     LocalView<int,2> matrix( array, shape );
/// which is identical for this matrix to previous Example 1
/// @endcode


template <typename Value, int Rank>
class LocalView {
    template <typename T>
    using is_non_const_value_type = typename std::is_same<T, typename std::remove_const<Value>::type>;

#define ENABLE_IF_NON_CONST                                                                             \
    template <bool EnableBool                                                                   = true, \
              typename std::enable_if<(!std::is_const<Value>::value && EnableBool), int>::type* = nullptr>

#define ENABLE_IF_CONST_WITH_NON_CONST(T)                                                                             \
    template <typename T, typename std::enable_if<(std::is_const<Value>::value && is_non_const_value_type<T>::value), \
                                                  int>::type* = nullptr>


public:
    // -- Type definitions
    using value_type  = Value;
    using return_type = value_type;

    static constexpr int RANK{Rank};

private:
    using slicer_t       = typename helpers::ArraySlicer<LocalView<Value, Rank>>;
    using const_slicer_t = typename helpers::ArraySlicer<const LocalView<const Value, Rank>>;

    template <typename... Args>
    struct slice_t {
        using type = typename slicer_t::template Slice<Args...>::type;
    };

    template <typename... Args>
    struct const_slice_t {
        using type = typename const_slicer_t::template Slice<Args...>::type;
    };

public:
    // -- Constructors


    template <typename ValueTp, typename = std::enable_if_t<std::is_convertible_v<ValueTp*, Value*>>>
    LocalView(const LocalView<ValueTp,Rank>& other): data_(other.data_), size_(other.size_), shape_(other.shape_), strides_(other.strides_) {}

    template <typename ValueTp, typename Int1, typename Int2, typename = std::enable_if_t<std::is_convertible_v<ValueTp*, Value*> && std::is_integral_v<Int1> && std::is_integral_v<Int2>>>
    LocalView(ValueTp* data, const Int1 shape[], const Int2 strides[]): data_(data) {
        size_ = 1;
        for (idx_t j = 0; j < Rank; ++j) {
            shape_[j]   = shape[j];
            strides_[j] = strides[j];
            size_ *= shape_[j];
        }
    }

    template <typename ValueTp, typename Int, typename = std::enable_if_t<std::is_convertible_v<ValueTp*, Value*> && std::is_integral_v<Int>>>
    LocalView(ValueTp* data, const Int shape[]): data_(data) {
        size_ = 1;
        for (int j = Rank - 1; j >= 0; --j) {
            shape_[j]   = shape[j];
            strides_[j] = size_;
            size_ *= shape_[j];
        }
    }

    template <typename ValueTp, typename ArrayShape, typename = std::enable_if_t<std::is_convertible_v<ValueTp*, Value*>>>
    LocalView(ValueTp* data, const ArrayShape& shape) : LocalView(data,shape.data()) {}

    ENABLE_IF_CONST_WITH_NON_CONST(value_type)
    operator const LocalView<value_type, Rank>&() const {
        static_assert(std::is_const<Value>::value, "must be const");
        static_assert(!std::is_const<value_type>::value, "must be non-const");
        return (const LocalView<value_type, Rank>&)(*this);
    }

    // -- Access methods

    template <typename... Idx, int Rank_ = Rank, typename = std::enable_if_t<sizeof...(Idx) == Rank_>>
    value_type& operator()(Idx... idx) {
        check_bounds(idx...);
        return data_[index(idx...)];
    }

    template <typename... Idx, int Rank_ = Rank, typename = std::enable_if_t<sizeof...(Idx) == Rank_>>
    const value_type& operator()(Idx... idx) const {
        check_bounds(idx...);
        return data_[index(idx...)];
    }

    template <typename Idx, int Rank_ = Rank, typename = std::enable_if_t<Rank_ == 1>>
    const value_type& operator[](Idx idx) const {
        check_bounds(idx);
        return data_[index(idx)];
    }

    template <typename Idx, int Rank_ = Rank, typename = std::enable_if_t<Rank_ == 1>>
    value_type& operator[](Idx idx) {
        check_bounds(idx);
        return data_[index(idx)];
    }

    idx_t size() const { return size_; }

    template <typename Int>
    idx_t shape(Int idx) const {
        return shape_[idx];
    }

    template <typename Int>
    idx_t stride(Int idx) const {
        return strides_[idx];
    }

    const idx_t* shape() const { return shape_.data(); }

    const idx_t* strides() const { return strides_.data(); }

    value_type const* data() const { return data_; }

    value_type* data() { return data_; }

    bool contiguous() const { return (size_ == shape_[0] * strides_[0] ? true : false); }

    ENABLE_IF_NON_CONST
    void assign(const value_type& value);

    void dump(std::ostream& os) const;

    static constexpr idx_t rank() { return Rank; }

    template <typename... Args>
    typename slice_t<Args...>::type slice(Args... args) {
        return slicer_t(*this).apply(args...);
    }

    template <typename... Args>
    typename const_slice_t<Args...>::type slice(Args... args) const {
        return const_slicer_t(*this).apply(args...);
    }

    friend std::ostream& operator<<(std::ostream& out, const LocalView& x) {
        x.dump(out);
        return out;
    }

private:
    // -- Private methods

    template <int Dim, typename Int, typename... Ints>
    constexpr idx_t index_part(Int idx, Ints... next_idx) const {
        return idx * strides_[Dim] + index_part<Dim + 1>(next_idx...);
    }

    template <int Dim, typename Int>
    constexpr idx_t index_part(Int last_idx) const {
        return last_idx * strides_[Dim];
    }

    template <typename... Ints>
    constexpr idx_t index(Ints... idx) const {
        return index_part<0>(idx...);
    }

#if ATLAS_ARRAYVIEW_BOUNDS_CHECKING
    template <typename... Ints>
    void check_bounds(Ints... idx) const {
        static_assert(sizeof...(idx) == Rank, "Expected number of indices is different from rank of array");
        return check_bounds_part<0>(idx...);
    }
#else
    template <typename... Ints>
    void check_bounds(Ints... idx) const {
        static_assert(sizeof...(idx) == Rank, "Expected number of indices is different from rank of array");
    }
#endif

    template <typename... Ints>
    void check_bounds_force(Ints... idx) const {
        static_assert(sizeof...(idx) == Rank, "Expected number of indices is different from rank of array");
        return check_bounds_part<0>(idx...);
    }

    template <int Dim, typename Int, typename... Ints>
    void check_bounds_part(Int idx, Ints... next_idx) const {
        if (idx_t(idx) >= shape_[Dim]) {
            throw_OutOfRange("LocalView", array_dim<Dim>(), idx, shape_[Dim]);
        }
        check_bounds_part<Dim + 1>(next_idx...);
    }

    template <int Dim, typename Int>
    void check_bounds_part(Int last_idx) const {
        if (idx_t(last_idx) >= shape_[Dim]) {
            throw_OutOfRange("LocalView", array_dim<Dim>(), last_idx, shape_[Dim]);
        }
    }

private:
    // -- Private data
    template<typename,int> friend class LocalView;

    value_type* data_;
    idx_t size_;
    std::array<idx_t,Rank> shape_;
    std::array<idx_t,Rank> strides_;

#undef ENABLE_IF_NON_CONST
#undef ENABLE_IF_CONST_WITH_NON_CONST
};

template<typename Value, int Rank>
using View = LocalView<Value,Rank>;

}  // namespace array
}  // namespace atlas
