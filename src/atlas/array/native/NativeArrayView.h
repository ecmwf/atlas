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

template <typename Value, int Rank, Intent AccessMode = Intent::ReadWrite>
class ArrayView {
public:
    // -- Type definitions
    using value_type = typename remove_const<Value>::type;
    using return_type =
        typename std::conditional<( AccessMode == Intent::ReadWrite ), value_type, value_type const>::type;

    static constexpr Intent ACCESS{AccessMode};
    static constexpr int RANK{Rank};

private:
    using slicer_t       = typename helpers::ArraySlicer<ArrayView<Value, Rank, AccessMode>>;
    using const_slicer_t = typename helpers::ArraySlicer<const ArrayView<Value, Rank, AccessMode>>;

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

    ArrayView( const ArrayView& other ) :
        data_( other.data_ ),
        size_( other.size_ ),
        shape_( other.shape_ ),
        strides_( other.strides_ ) {}

    ArrayView( const value_type* data, const ArrayShape& shape, const ArrayStrides& strides ) :
        data_( const_cast<value_type*>( data ) ) {
        size_ = 1;
        for ( int j = 0; j < Rank; ++j ) {
            shape_[j]   = shape[j];
            strides_[j] = strides[j];
            size_ *= shape_[j];
        }
    }

    // -- Access methods

    template <typename... Idx>
    return_type& operator()( Idx... idx ) {
        check_bounds( idx... );
        return data_[index( idx... )];
    }

    template <typename... Ints>
    const value_type& operator()( Ints... idx ) const {
        return data_[index( idx... )];
    }

    template <typename Int, bool EnableBool = true>
    typename std::enable_if<( Rank == 1 && EnableBool ), const value_type&>::type operator[]( Int idx ) const {
        check_bounds( idx );
        return data_[idx * strides_[0]];
    }

    template <typename Int, bool EnableBool = true>
    typename std::enable_if<( Rank == 1 && EnableBool ), value_type&>::type operator[]( Int idx ) {
        check_bounds( idx );
        return data_[idx * strides_[0]];
    }

    template <unsigned int Dim>
    idx_t shape() const {
        return shape( Dim );
    }

    template <unsigned int Dim>
    idx_t stride() const {
        return stride( Dim );
    }

    idx_t size() const { return size_; }

    static constexpr idx_t rank() { return Rank; }

    const idx_t* strides() const { return strides_.data(); }

    const idx_t* shape() const { return shape_.data(); }

    template <typename Int>
    idx_t shape( Int idx ) const {
        return shape_[idx];
    }

    template <typename Int>
    idx_t stride( Int idx ) const {
        return strides_[idx];
    }

    value_type const* data() const { return data_; }

    value_type* data() { return data_; }

    bool valid() const { return true; }

    bool contiguous() const { return ( size_ == shape_[0] * strides_[0] ? true : false ); }

    void assign( const value_type& value );

    void assign( const std::initializer_list<value_type>& list );

    void dump( std::ostream& os ) const;

    template <typename... Args>
    typename slice_t<Args...>::type slice( Args... args ) {
        return slicer_t( *this ).apply( args... );
    }

    template <typename... Args>
    typename const_slice_t<Args...>::type slice( Args... args ) const {
        return const_slicer_t( *this ).apply( args... );
    }

private:
    // -- Private methods

    template <int Dim, typename Int, typename... Ints>
    constexpr idx_t index_part( Int idx, Ints... next_idx ) const {
        return idx * strides_[Dim] + index_part<Dim + 1>( next_idx... );
    }

    template <int Dim, typename Int>
    constexpr idx_t index_part( Int last_idx ) const {
        return last_idx * strides_[Dim];
    }

    template <typename... Ints>
    constexpr idx_t index( Ints... idx ) const {
        return index_part<0>( idx... );
    }

#if ATLAS_ARRAYVIEW_BOUNDS_CHECKING
    template <typename... Ints>
    void check_bounds( Ints... idx ) const {
        static_assert( sizeof...( idx ) == Rank, "Expected number of indices is different from rank of array" );
        return check_bounds_part<0>( idx... );
    }
#else
    template <typename... Ints>
    void check_bounds( Ints... idx ) const {
        static_assert( sizeof...( idx ) == Rank, "Expected number of indices is different from rank of array" );
    }
#endif

    template <typename... Ints>
    void check_bounds_force( Ints... idx ) const {
        static_assert( sizeof...( idx ) == Rank, "Expected number of indices is different from rank of array" );
        return check_bounds_part<0>( idx... );
    }

    template <int Dim, typename Int, typename... Ints>
    void check_bounds_part( Int idx, Ints... next_idx ) const {
        if ( idx_t( idx ) >= shape_[Dim] ) {
            throw_OutOfRange( "ArrayView", array_dim<Dim>(), idx, shape_[Dim] );
        }
        check_bounds_part<Dim + 1>( next_idx... );
    }

    template <int Dim, typename Int>
    void check_bounds_part( Int last_idx ) const {
        if ( idx_t( last_idx ) >= shape_[Dim] ) {
            throw_OutOfRange( "ArrayView", array_dim<Dim>(), last_idx, shape_[Dim] );
        }
    }

    // -- Private data

    value_type* data_;
    idx_t size_;
    std::array<idx_t, Rank> shape_;
    std::array<idx_t, Rank> strides_;
};

//------------------------------------------------------------------------------------------------------

}  // namespace array
}  // namespace atlas
