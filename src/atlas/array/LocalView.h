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
#include "atlas/mdspan.h"

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


template <typename ElementType, int Rank>
class LocalView {
    template <typename T>
    using is_non_const_value_type = typename std::is_same<T, typename std::remove_const<ElementType>::type>;

#define ENABLE_IF_NON_CONST                                                                             \
    template <bool EnableBool                                                                   = true, \
              typename std::enable_if<(!std::is_const<ElementType>::value && EnableBool), int>::type* = nullptr>

#define ENABLE_IF_CONST_WITH_NON_CONST(T)                                                                             \
    template <typename T, typename std::enable_if<(std::is_const<ElementType>::value && is_non_const_value_type<T>::value), \
                                                  int>::type* = nullptr>


public:
    // -- Type definitions
    using element_type = ElementType;
    using value_type  = std::remove_cv_t<element_type>;
    using return_type = element_type;
    using data_handle_type = element_type*;

    static constexpr int RANK{Rank};

private:
    using slicer_t       = typename helpers::ArraySlicer<LocalView<ElementType, Rank>>;
    using const_slicer_t = typename helpers::ArraySlicer<const LocalView<const ElementType, Rank>>;

    template <typename... Args>
    struct slice_t {
        using type = typename slicer_t::template Slice<Args...>::type;
    };

    template <typename... Args>
    struct const_slice_t {
        using type = typename const_slicer_t::template Slice<Args...>::type;
    };

    using mdspan_extents_type = dextents<size_t,RANK>;
    using mdspan_strides_type = std::array<size_t,RANK>;

public:
    // -- Constructors

    LocalView(const LocalView& other): data_(other.data_), size_(other.size_) {
        init_metadata_pointers();
        copy_metadata(other.shape_, other.strides_);
    }


    template <typename ElementTypeTp, typename = std::enable_if_t<std::is_convertible_v<ElementTypeTp*, ElementType*>>>
    LocalView(const LocalView<ElementTypeTp,Rank>& other): data_(other.data_), size_(other.size_) {
        init_metadata_pointers();
        copy_metadata(other.shape_, other.strides_);
    }

    LocalView& operator=(const LocalView& other) {
        if (this != &other) {
            data_ = other.data_;
            size_ = other.size_;
            init_metadata_pointers();
            copy_metadata(other.shape_, other.strides_);
        }
        return *this;
    }

    template <typename ElementTypeTp, typename = std::enable_if_t<std::is_convertible_v<ElementTypeTp*, ElementType*>>>
    LocalView& operator=(const LocalView<ElementTypeTp,Rank>& other) {
        data_ = other.data_;
        size_ = other.size_;
        init_metadata_pointers();
        copy_metadata(other.shape_, other.strides_);
        return *this;
    }

    template <typename ElementTypeTp, typename Int1, typename Int2, typename = std::enable_if_t<std::is_convertible_v<ElementTypeTp*, ElementType*> && std::is_integral_v<Int1> && std::is_integral_v<Int2>>>
    LocalView(ElementTypeTp* data, const Int1 shape[], const Int2 strides[]): data_(data) {
        init_metadata_pointers();
        size_ = 1;
        for (idx_t j = 0; j < Rank; ++j) {
            shape_[j]   = shape[j];
            strides_[j] = strides[j];
            size_ *= shape_[j];
        }
    }

    template <typename ElementTypeTp, typename Int, typename = std::enable_if_t<std::is_convertible_v<ElementTypeTp*, ElementType*> && std::is_integral_v<Int>>>
    LocalView(ElementTypeTp* data, const Int shape[]): data_(data) {
        init_metadata_pointers();
        size_ = 1;
        for (int j = Rank - 1; j >= 0; --j) {
            shape_[j]   = shape[j];
            strides_[j] = size_;
            size_ *= shape_[j];
        }
    }

    template <typename ElementTypeTp, typename ArrayShape, typename = std::enable_if_t<std::is_convertible_v<ElementTypeTp*, ElementType*>>>
    LocalView(ElementTypeTp* data, const ArrayShape& shape) : LocalView(data,shape.data()) {}


    template <typename T, typename E, typename L, typename A, typename = std::enable_if_t<std::is_convertible_v<typename A::data_handle_type, ElementType*> && E::rank() == Rank>>
    LocalView(mdspan<T,E,L,A>& other) :
        data_(other.data_handle()), size_(other.size()) {
        init_metadata_pointers();
        for (int j = 0; j < Rank; ++j) {
            shape_[j] = other.extent(j);
            strides_[j] = other.stride(j);
        }
    }

    ENABLE_IF_CONST_WITH_NON_CONST(element_type)
    operator const LocalView<element_type, Rank>&() const {
        static_assert(std::is_const<element_type>::value, "must be const");
        static_assert(!std::is_const<value_type>::value, "must be non-const");
        return (const LocalView<element_type, Rank>&)(*this);
    }


    // -- Access methods

    template <typename... Idx, int Rank_ = Rank, typename = std::enable_if_t<sizeof...(Idx) == Rank_>>
    inline ATLAS_HOST_DEVICE
    element_type& operator()(Idx... idx) {
        check_bounds(idx...);
        return data_[index(idx...)];
    }

    template <typename... Idx, int Rank_ = Rank, typename = std::enable_if_t<sizeof...(Idx) == Rank_>>
    inline ATLAS_HOST_DEVICE
    const element_type& operator()(Idx... idx) const {
        check_bounds(idx...);
        return data_[index(idx...)];
    }

    template <typename Idx, int Rank_ = Rank, typename = std::enable_if_t<Rank_ == 1>>
    inline ATLAS_HOST_DEVICE
    const element_type& operator[](Idx idx) const noexcept(!ATLAS_ARRAYVIEW_BOUNDS_CHECKING || !ATLAS_HOST_COMPILE) {
        check_bounds(idx);
        return data_[index(idx)];
    }

    template <typename Idx, int Rank_ = Rank, typename = std::enable_if_t<Rank_ == 1>>
    inline ATLAS_HOST_DEVICE
    element_type& operator[](Idx idx) noexcept(!ATLAS_ARRAYVIEW_BOUNDS_CHECKING || !ATLAS_HOST_COMPILE) {
        check_bounds(idx);
        return data_[index(idx)];
    }

    inline ATLAS_HOST_DEVICE
    idx_t size() const noexcept { return size_; }

    template <typename Int>
    inline ATLAS_HOST_DEVICE
    idx_t shape(Int idx) const noexcept {
        return shape_[idx];
    }

    /// @brief Return number of values in dimension idx, equivalent to shape(idx)
    template <typename Int>
    inline ATLAS_HOST_DEVICE
    idx_t extent(Int idx) const noexcept {
        return shape(idx);
    }

    template <typename Int>
    inline ATLAS_HOST_DEVICE
    idx_t stride(Int idx) const noexcept {
        return strides_[idx];
    }

    inline ATLAS_HOST_DEVICE
    const idx_t* shape() const noexcept { return shape_; }

    inline ATLAS_HOST_DEVICE
    const idx_t* strides() const noexcept { return strides_; }

    inline ATLAS_HOST_DEVICE
    element_type const* data() const noexcept { return data_; }

    inline ATLAS_HOST_DEVICE
    element_type* data() noexcept { return data_; }

    inline ATLAS_HOST_DEVICE
    constexpr data_handle_type data_handle() const noexcept { return data_; }

    inline ATLAS_HOST_DEVICE
    bool contiguous() const noexcept { return (size_ == shape_[0] * strides_[0] ? true : false); }

    ENABLE_IF_NON_CONST
    void assign(const value_type& value);

    void dump(std::ostream& os) const;

    inline static constexpr idx_t rank() noexcept { return Rank; }

    template <typename... Args>
    inline typename slice_t<Args...>::type slice(Args... args) {
        return slicer_t(*this).apply(args...);
    }

    template <typename... Args>
    inline typename const_slice_t<Args...>::type slice(Args... args) const {
        return const_slicer_t(*this).apply(args...);
    }

    friend std::ostream& operator<<(std::ostream& out, const LocalView& x) {
        x.dump(out);
        return out;
    }

private:

    inline ATLAS_HOST_DEVICE 
    void init_metadata_pointers() noexcept {
        shape_ = shape_data_;
        strides_ = strides_data_;
    }

    template <typename Shape, typename Strides>
    inline ATLAS_HOST_DEVICE
    void copy_metadata(const Shape& shape, const Strides& strides) noexcept {
        for (int j = 0; j < Rank; ++j) {
            shape_[j] = shape[j];
            strides_[j] = strides[j];
        }
    }

    // -- Private methods

    template <int Dim, typename Int, typename... Ints>
    inline ATLAS_HOST_DEVICE
    constexpr idx_t index_part(Int idx, Ints... next_idx) const {
        return idx * strides_[Dim] + index_part<Dim + 1>(next_idx...);
    }

    template <int Dim, typename Int>
    inline ATLAS_HOST_DEVICE
    constexpr idx_t index_part(Int last_idx) const {
        return last_idx * strides_[Dim];
    }

    template <typename... Ints>
    inline ATLAS_HOST_DEVICE
    constexpr idx_t index(Ints... idx) const {
        return index_part<0>(idx...);
    }

#if ATLAS_ARRAYVIEW_BOUNDS_CHECKING
    template <typename... Ints>
    inline ATLAS_HOST_DEVICE
    void check_bounds(Ints... idx) const {
        static_assert(sizeof...(idx) == Rank, "Expected number of indices is different from rank of array");
#if ATLAS_HOST_COMPILE
        return check_bounds_part<0>(idx...);
#endif
    }
#else
    template <typename... Ints>
    inline ATLAS_HOST_DEVICE
    void check_bounds(Ints... idx) const noexcept {
        static_assert(sizeof...(idx) == Rank, "Expected number of indices is different from rank of array");
    }
#endif

    template <typename... Ints>
    inline ATLAS_HOST_DEVICE
    void check_bounds_force(Ints... idx) const {
        static_assert(sizeof...(idx) == Rank, "Expected number of indices is different from rank of array");
#if ATLAS_HOST_COMPILE
        return check_bounds_part<0>(idx...);
#endif
    }

    template <int Dim, typename Int, typename... Ints>
    inline void check_bounds_part(Int idx, Ints... next_idx) const {
        if (idx_t(idx) >= shape_[Dim]) {
            throw_OutOfRange("LocalView", array_dim<Dim>(), idx, shape_[Dim]);
        }
        check_bounds_part<Dim + 1>(next_idx...);
    }

    template <int Dim, typename Int>
    inline void check_bounds_part(Int last_idx) const {
        if (idx_t(last_idx) >= shape_[Dim]) {
            throw_OutOfRange("LocalView", array_dim<Dim>(), last_idx, shape_[Dim]);
        }
    }

    template<int... i>
    inline ATLAS_HOST_DEVICE
    mdspan_extents_type _get_mdspan_extents(std::integer_sequence<int, i...> = {}) const {
        return mdspan_extents_type{(shape_[i])...};
    }

    inline ATLAS_HOST_DEVICE
    mdspan_extents_type mdspan_extents() const {
        return _get_mdspan_extents(std::make_integer_sequence<int,RANK>{});
    }

    template<int... i>
    inline ATLAS_HOST_DEVICE
    mdspan_strides_type _get_mdspan_strides(std::integer_sequence<int, i...> = {}) const {
        return mdspan_strides_type{(static_cast<typename mdspan_strides_type::value_type>(strides_[i]))...};
    }

    inline ATLAS_HOST_DEVICE
    mdspan_strides_type mdspan_strides() const {
        return _get_mdspan_strides(std::make_integer_sequence<int,RANK>{});
    }

private:
    // -- Private data
    template<typename,int> friend class LocalView;

    element_type* data_;
    idx_t size_;
    idx_t* shape_;
    idx_t* strides_;
    idx_t shape_data_[Rank];
    idx_t strides_data_[Rank];

#undef ENABLE_IF_NON_CONST
#undef ENABLE_IF_CONST_WITH_NON_CONST
};

template<typename Value, int Rank>
using View = LocalView<Value,Rank>;

}  // namespace array
}  // namespace atlas
