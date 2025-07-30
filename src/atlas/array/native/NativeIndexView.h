/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

/// @file IndexView.h
/// This file contains the IndexView class, a class that allows to wrap any
/// contiguous raw data into
/// a view which is accessible with multiple indices.
/// This view is intended to work with Connectivity Tables storing Fortran
/// Numbering internally
/// All it needs is the strides for each index, and the shape of each index.
/// ATTENTION: The last index is stride 1
///
/// Bounds-checking can be turned ON by defining
/// "ATLAS_INDEXVIEW_BOUNDS_CHECKING"
/// before including this header.
///
/// Example:
/// int[] array = { 1, 2, 3, 4, 5, 6, 7, 8, 9};
/// int[2] strides = { 3, 1 };
/// int[2] shape = { 3, 3 };
/// IndexView<int,2> matrix( array, strides, shape );
/// for( size_t i=0; i<matrix.shape(0); ++i ) {
///   for( size_t j=0; j<matrix.shape(1); ++j ) {
///     matrix(i,j) *= 10;
///   }
/// }
///
/// There is also an easier way to wrap Field and Array classes:
/// IndexView<int,3> fieldview( Field );
/// IndexView<int,2> INDEXVIEW( Array );

#pragma once

#include <initializer_list>
#include <iosfwd>
#include <type_traits>
#include <utility>

#include "atlas/array/ArrayDataStore.h"
#include "atlas/library/config.h"
#include "atlas/mdspan.h"

//------------------------------------------------------------------------------------------------------

namespace atlas {
namespace array {

//------------------------------------------------------------------------------------------------------

#if ATLAS_HAVE_FORTRAN
#define INDEX_REF Index
#define FROM_FORTRAN -1
#else
#define INDEX_REF *
#define FROM_FORTRAN
#endif

#define ENABLE_IF_NON_CONST                                                                             \
    template <bool EnableBool                                                                   = true, \
              typename std::enable_if<(!std::is_const<Value>::value && EnableBool), int>::type* = nullptr>

//------------------------------------------------------------------------------------------------------

///@brief Multidimensional access to Array or Field objects containing index fields
///       that are compatible with Fortran indexing
///
/// Index fields that are compatible with Fortran are 1-based: a value `1` corresponds
/// to the first index in a Fortran array. In C++ the first index would of course correspond to
/// the value `0`. To stay compatible, we created this class IndexView that subtracts `1`
/// when an index value is accessed, and adds `1` when a value is set.
///
/// If Atlas is compiled without the FORTRAN feature (ATLAS_HAVE_FORTRAN==0), then the addition
/// and substraction of `1` is compiled out, which may slightly improve performance.
template <typename Value, int Rank>
class IndexView {
public:
    using value_type = typename remove_const<Value>::type;
    using accessor_type = index_accessor<Value,ATLAS_HAVE_FORTRAN>;
    using reference = typename accessor_type::reference;

private:
    using mdspan_extents_type = dextents<size_t,Rank>;
    using mdspan_strides_type = std::array<size_t,Rank>;
    using mdspan_type         = mdspan<Value, mdspan_extents_type, layout_stride, index_accessor<Value,ATLAS_HAVE_FORTRAN>>;
    using const_mdspan_type   = mdspan<const Value, mdspan_extents_type, layout_stride, index_accessor<const Value,ATLAS_HAVE_FORTRAN>>;

public:
#ifndef DOXYGEN_SHOULD_SKIP_THIS
    IndexView(Value* data, const idx_t shape[Rank]);

    IndexView(Value* data, const idx_t shape[Rank], const idx_t strides[Rank]);
#endif

    template <typename T, typename E, typename L, typename A, typename = std::enable_if_t<std::is_convertible_v<typename A::data_handle_type, Value*> && E::rank() == Rank>>
    IndexView(mdspan<T,E,L,A>& other) :
        data_(other.data_handle()) {
        for (int j = 0; j < Rank; ++j) {
            shape_[j] = other.extent(j);
            strides_[j] = other.stride(j);
        }
    }

    // -- Access methods

    /// @brief Multidimensional index operator: view(i,j,k,...)
    template <typename... Idx>
    reference operator()(Idx... idx) {
        check_bounds(idx...);
        return accessor_.access(data_, index(idx...));
    }

    template <typename... Idx>
    value_type operator()(Idx... idx) const {
        check_bounds(idx...);
        return data_[index(idx...)] FROM_FORTRAN;
    }

    ENABLE_IF_NON_CONST
    void assign(const std::initializer_list<value_type>& list);

    template <typename Int>
    idx_t shape(Int idx) const {
        return shape_[idx];
    }

    mdspan_type as_mdspan() {
        return mdspan_type{data_, {mdspan_extents(), mdspan_strides()}};
    }

    const_mdspan_type as_mdspan() const {
        return const_mdspan_type{data_, {mdspan_extents(), mdspan_strides()}};
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

#if ATLAS_INDEXVIEW_BOUNDS_CHECKING
    template <typename... Ints>
    void check_bounds(Ints... idx) const {
        static_assert(sizeof...(idx) == Rank, "Expected number of indices is different from rank of array");
        return check_bounds_part<0>(idx...);
    }
#else
    template <typename... Ints>
    void check_bounds(Ints...) const {}
#endif

    template <typename... Ints>
    void check_bounds_force(Ints... idx) const {
        static_assert(sizeof...(idx) == Rank, "Expected number of indices is different from rank of array");
        return check_bounds_part<0>(idx...);
    }

    template <int Dim, typename Int, typename... Ints>
    void check_bounds_part(Int idx, Ints... next_idx) const {
        if (idx_t(idx) >= shape_[Dim]) {
            throw_OutOfRange("IndexView", array_dim<Dim>(), idx, shape_[Dim]);
        }
        check_bounds_part<Dim + 1>(next_idx...);
    }

    template <int Dim, typename Int>
    void check_bounds_part(Int last_idx) const {
        if (idx_t(last_idx) >= shape_[Dim]) {
            throw_OutOfRange("IndexView", array_dim<Dim>(), last_idx, shape_[Dim]);
        }
    }

    idx_t size() const { return shape_[0]; }

    void dump(std::ostream& os) const;

    template<int... i>
    mdspan_extents_type _get_mdspan_extents(std::integer_sequence<int, i...> = {}) const {
        return mdspan_extents_type{(shape_[i])...};
    }

    mdspan_extents_type mdspan_extents() const {
        return _get_mdspan_extents(std::make_integer_sequence<int,Rank>{});
    }

    template<int... i>
    mdspan_strides_type _get_mdspan_strides(std::integer_sequence<int, i...> = {}) const {
        return mdspan_strides_type{(static_cast<typename mdspan_strides_type::value_type>(strides_[i]))...};
    }

    mdspan_strides_type mdspan_strides() const {
        return _get_mdspan_strides(std::make_integer_sequence<int,Rank>{});
    }

private:
    Value* data_;
    idx_t strides_[Rank];
    idx_t shape_[Rank];
    static constexpr accessor_type accessor_{};
};

template <typename Value, int Rank>
class LocalIndexView : public IndexView<Value, Rank> {
    using Base = IndexView<Value, Rank>;

public:
    using Base::Base;
};

#undef INDEX_REF
#undef FROM_FORTRAN
#undef ENABLE_IF_NON_CONST

//------------------------------------------------------------------------------------------------------

}  // namespace array
}  // namespace atlas
