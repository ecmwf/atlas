/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

/// @file ArrayView.h
/// This file contains the ArrayView class, a class that allows to wrap any
/// contiguous raw data into
/// a view which is accessible with multiple indices.
/// All it needs is the strides for each index, and the shape of each index.
/// ATTENTION: The last index is stride 1
///
/// Bounds-checking can be turned ON by defining
/// "ATLAS_ARRAYVIEW_BOUNDS_CHECKING"
/// before including this header.
///
/// Example 1:
///     int[] array = { 1, 2, 3, 4, 5, 6, 7, 8, 9};
///     int[2] strides = { 3, 1 };
///     int[2] shape = { 3, 3 };
///     ArrayView<int,2> matrix( array, shape, strides );
///     for( idx_t i=0; i<matrix.shape(0); ++i ) {
///       for( idx_t j=0; j<matrix.shape(1); ++j ) {
///         matrix(i,j) *= 10;
///       }
///     }
///
/// Strides can also be omitted as for most common cases it can be inferred
/// from the shape.
///
/// Example 2:
///     int[] array = { 1, 2, 3, 4, 5, 6, 7, 8, 9};
///     int[2] shape = { 3, 3 };
///     ArrayView<int,2> matrix( array, shape );
/// which is identical for this matrix to previous Example 1
///
/// There is also an easier way to wrap Field and Array classes:
///
/// Example 3:
///     ArrayView<int,3> fieldview( Field );
///     ArrayView<int,2> arrayview( Array );
///
/// @author Willem Deconinck

#pragma once

#include <array>
#include <cstddef>
#include <initializer_list>
#include <iostream>
#include <type_traits>

#include "atlas/array/ArrayUtil.h"
#include "atlas/array/ArrayViewDefs.h"
#include "atlas/array/LocalView.h"
#include "atlas/array/Range.h"
#include "atlas/array/helpers/ArraySlicer.h"
#include "atlas/library/config.h"

//------------------------------------------------------------------------------------------------------

namespace atlas {
namespace array {

//------------------------------------------------------------------------------------------------------

/// @brief Multi-dimensional access to a Array object or Field object
///
/// An ArrayView enables access to the inner data-memory-storage of the Array.
/// It is required to create the view using the make_view helper functions.
///
/// ### Example 1:
///
/// @code{.cpp}
///    auto view = make_view<double,2>( Array );
///    double sum = 0;
///    for( idx_t i=0; i<view.shape(0); ++i ) {
///       for( idx_t j=0; j<view.shape(1); ++j ) {
///           sum += view(i,j);
///       }
///    }
/// @endcode
///
/// ### Data storage
///
/// Depending on whether atlas was compiled with the Feature `GRIDTOOLS_STORAGE`, the
/// internal data allocation of the Array object is managed by [GridTools](https://gridtools.github.io/)
/// (ATLAS_HAVE_GRIDTOOLS==1) or by Atlas itself (ATLAS_HAVE_GRIDTOOLS==0).
///
/// The ArrayView class is therefore also compiled differently dependening on this feature.

template <typename Value, int Rank>
class ArrayView {
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
    using value_type                   = Value;
    using non_const_value_type         = typename std::remove_const<Value>::type;
    static constexpr bool is_const     = std::is_const<Value>::value;
    static constexpr bool is_non_const = !std::is_const<Value>::value;
    static constexpr int RANK{Rank};

private:
    using slicer_t       = typename helpers::ArraySlicer<ArrayView<Value, Rank>>;
    using const_slicer_t = typename helpers::ArraySlicer<const ArrayView<const Value, Rank>>;

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

    ArrayView(const ArrayView& other):
        data_(other.data_), size_(other.size_), shape_(other.shape_), strides_(other.strides_) {}

    ENABLE_IF_CONST_WITH_NON_CONST(value_type)
    ArrayView(const ArrayView<value_type, Rank>& other): data_(other.data()), size_(other.size()) {
        for (idx_t j = 0; j < Rank; ++j) {
            shape_[j]   = other.shape(j);
            strides_[j] = other.stride(j);
        }
    }

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    // This constructor should not be used directly, but only through a array::make_view() function.
    ArrayView(value_type* data, const ArrayShape& shape, const ArrayStrides& strides): data_(data) {
        size_ = 1;
        for (int j = 0; j < Rank; ++j) {
            shape_[j]   = shape[j];
            strides_[j] = strides[j];
            size_ *= size_t(shape_[j]);
        }
    }
#endif

    ENABLE_IF_CONST_WITH_NON_CONST(value_type)
    operator const ArrayView<value_type, Rank>&() const { return *(const ArrayView<value_type, Rank>*)(this); }


    // -- Access methods

    /// @brief Multidimensional index operator: view(i,j,k,...)
    template <typename... Idx>
    value_type& operator()(Idx... idx) {
        check_bounds(idx...);
        return data_[index(idx...)];
    }

    /// @brief Multidimensional index operator: view(i,j,k,...)
    template <typename... Ints>
    const value_type& operator()(Ints... idx) const {
        return data_[index(idx...)];
    }

    /// @brief Access to data using square bracket [idx] operator @m_class{m-label m-warning} **Rank==1**.
    ///
    /// Note that this function is only present when Rank == 1
#ifndef DOXYGEN_SHOULD_SKIP_THIS
    template <typename Int, bool EnableBool = true>
    typename std::enable_if<(Rank == 1 && EnableBool), const value_type&>::type operator[](Int idx) const {
#else
    // Doxygen API is cleaner!
    template <typename Int>
    value_type operator[](Int idx) const {
#endif
        check_bounds(idx);
        return data_[idx * strides_[0]];
    }

    /// @brief Access to data using square bracket [idx] operator @m_class{m-label m-warning} **Rank==1**
    ///
    /// Note that this function is only present when Rank == 1
#ifndef DOXYGEN_SHOULD_SKIP_THIS
    template <typename Int, bool EnableBool = true>
    typename std::enable_if<(Rank == 1 && EnableBool), value_type&>::type operator[](Int idx) {
#else
    // Doxygen API is cleaner!
    template <typename Int>
    value_type operator[](Int idx) {
#endif
        check_bounds(idx);
        return data_[idx * strides_[0]];
    }

    /// @brief Return number of values in dimension **Dim** (template argument)
    ///
    /// Example use:
    /// @code{.cpp}
    /// for( idx_t i=0; i<view.shape<0>(); ++i ) {
    ///   for( idx_t j=0; j<view.shape<1>(); +j ) {
    ///     ...
    ///   }
    /// }
    /// @endcode
    template <unsigned int Dim>
    idx_t shape() const {
        return shape_[Dim];
    }

    /// @brief Return stride for values in dimension **Dim** (template argument)
    template <unsigned int Dim>
    idx_t stride() const {
        return strides_[Dim];
    }

    /// @brief Return total number of values (accumulated over all dimensions)
    size_t size() const { return size_; }

    /// @brief Return the number of dimensions
    static constexpr idx_t rank() { return Rank; }

    const idx_t* strides() const { return strides_.data(); }

    const idx_t* shape() const { return shape_.data(); }

    /// @brief Return number of values in dimension idx
    template <typename Int>
    idx_t shape(Int idx) const {
        return shape_[idx];
    }

    /// @brief Return stride for values in dimension idx
    template <typename Int>
    idx_t stride(Int idx) const {
        return strides_[idx];
    }

    /// @brief Access to internal data. @m_class{m-label m-danger} **dangerous**
    value_type const* data() const { return data_; }

    /// @brief Access to internal data. @m_class{m-label m-danger} **dangerous**
    value_type* data() { return data_; }

    bool valid() const { return true; }

    /// @brief Return true when all values are contiguous in memory.
    ///
    /// This means that if there is e.g. padding in the fastest dimension, or if
    /// the ArrayView represents a slice, the returned value will be false.
    bool contiguous() const { return (size_ == size_t(shape_[0]) * size_t(strides_[0]) ? true : false); }

    ENABLE_IF_NON_CONST
    void assign(const value_type& value);

    ENABLE_IF_NON_CONST
    void assign(const std::initializer_list<value_type>& list);

    ENABLE_IF_NON_CONST
    void assign(const ArrayView& other);

    void dump(std::ostream& os) const;

    /// @brief Obtain a slice from this view:  view.slice( Range, Range, ... )
    ///
    /// The return type of this function is intentionally `auto` and is guaranteed to have
    /// the same API as ArrayView, but is not necessarily this type.
    ///
    /// If the current view has Rank == 2, a Rank == 1 slice can be created in several ways:
    ///
    /// @code{.cpp}
    ///   auto slice1 = view.slice( Range(0,2), Range::all() );
    ///   auto slice2 = view.slice( Range::To(2), Range::all() );
    /// @endcode
    ///
    /// Sometimes it may be required to extend the rank of the current view to cater for
    /// certain algorithms requiring an extra rank.
    ///
    /// @code{.cpp}
    ///   auto slice3 = view.slice( Range::all(), Range::all(), Range::dummy() );
    /// @endcode
    template <typename... Args>
    auto slice(Args... args) {
        return slicer_t(*this).apply(args...);
    }


    /// @brief Obtain a slice from this view:  view.slice( Range, Range, ... )
    template <typename... Args>
    auto slice(Args... args) const {
        return const_slicer_t(*this).apply(args...);
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
            throw_OutOfRange("ArrayView", array_dim<Dim>(), idx, shape_[Dim]);
        }
        check_bounds_part<Dim + 1>(next_idx...);
    }

    template <int Dim, typename Int>
    void check_bounds_part(Int last_idx) const {
        if (idx_t(last_idx) >= shape_[Dim]) {
            throw_OutOfRange("ArrayView", array_dim<Dim>(), last_idx, shape_[Dim]);
        }
    }

    // -- Private data

    value_type* data_;
    size_t size_;
    std::array<idx_t, Rank> shape_;
    std::array<idx_t, Rank> strides_;
};

//------------------------------------------------------------------------------------------------------

#undef ENABLE_IF_NON_CONST
#undef ENABLE_IF_CONST_WITH_NON_CONST

}  // namespace array
}  // namespace atlas
