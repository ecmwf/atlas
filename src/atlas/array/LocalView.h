/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

/// @file LocalView.h
/// This file contains the LocalView class, a class that allows to wrap any contiguous raw data into
/// a view which is accessible with multiple indices.
/// All it needs is the strides for each index, and the shape of each index.
/// ATTENTION: The last index is stride 1
///
/// Bounds-checking can be turned ON by defining "ATLAS_ARRAYVIEW_BOUNDS_CHECKING"
/// before including this header.
///
/// Example 1:
///     int[] array = { 1, 2, 3, 4, 5, 6, 7, 8, 9};
///     int[2] strides = { 3, 1 };
///     int[2] shape = { 3, 3 };
///     LocalView<int,2> matrix( array, shape, strides );
///     for( size_t i=0; i<matrix.shape(0); ++i ) {
///       for( size_t j=0; j<matrix.shape(1); ++j ) {
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
///     LocalView<int,2> matrix( array, shape );
/// which is identical for this matrix to previous Example 1
///
/// @author Willem Deconinck
#pragma once

#include <cstddef>
#include <type_traits>
#include "atlas/library/config.h"
#include "atlas/array/ArrayUtil.h"

//------------------------------------------------------------------------------------------------------

namespace atlas {
namespace array {

template< typename Value, int Rank >
class LocalView {
public:

// -- Type definitions
    using value_type = typename remove_const<Value>::type;
    using Slice = typename std::conditional<(Rank==1), value_type&, LocalView<value_type,Rank-1> >::type;

public:

// -- Constructors

    LocalView( const value_type* data, const size_t shape[], const size_t strides[] ) :
        data_(const_cast<value_type*>(data)) {
        size_ = 1;
        for( size_t j=0; j<Rank; ++j ) {
            shape_[j] = shape[j];
            strides_[j] = strides[j];
            size_ *= shape_[j];
        }
    }

    LocalView( const value_type* data, const size_t shape[] ) :
        data_(const_cast<value_type*>(data)) {
        size_ = 1;
        for( int j=Rank-1; j>=0; --j ) {
            shape_[j] = shape[j];
            strides_[j] = size_;
            size_ *= shape_[j];
        }
    }

    LocalView( const value_type* data, const ArrayShape& shape ) :
        data_(const_cast<value_type*>(data)) {
        size_ = 1;
        for( int j=Rank-1; j>=0; --j ) {
            shape_[j]   = shape[j];
            strides_[j] = size_;
            size_ *= shape_[j];
        }
    }

// -- Access methods

    template < typename... Ints >
    value_type& operator()(Ints... idx) {
        check_bounds(idx...);
        return data_[index(idx...)];
    }

    template < typename... Ints >
    const value_type& operator()(Ints... idx) const {
        check_bounds(idx...);
        return data_[index(idx...)];
    }

    Slice at(const size_t i) const {
        return Slicer<Slice, Rank==1>(*this).apply(i);
    }

    size_t size() const { return size_;}

    size_t shape(size_t idx) const { return shape_[idx]; }

    value_type const* data() const { return data_; }

    value_type*       data()       { return data_; }

    bool contiguous() const {
        return (size_ == shape_[0]*strides_[0] ? true : false);
    }

    void assign(const value_type& value);

    void dump(std::ostream& os) const;

    static constexpr int rank() { return Rank; }

private:

// -- Type definitions

    template <typename ReturnType = Slice, bool ToScalar = false>
    struct Slicer {
        Slicer(LocalView<value_type, Rank> const& lv) : lv_(lv) {}
        LocalView<value_type, Rank> const& lv_;
        ReturnType apply(const size_t i) const {
            return LocalView<value_type, Rank - 1>(lv_.data_ + lv_.strides_[0] * i, lv_.shape_ + 1, lv_.strides_ + 1);
        }
    };

    template <typename ReturnType>
    struct Slicer<ReturnType, true> {
        Slicer(LocalView<value_type, Rank> const& lv) : lv_(lv) {}
        LocalView<value_type, Rank> const& lv_;
        ReturnType apply(const size_t i) const {
            return *(lv_.data_ + lv_.strides_[0] * i);
        }
    };

// -- Private methods

    template < int Dim, typename Int, typename... Ints >
    constexpr int index_part(Int idx, Ints... next_idx) const {
        return idx*strides_[Dim] + index_part<Dim+1>( next_idx... );
    }

    template < int Dim, typename Int >
    constexpr int index_part(Int last_idx) const {
        return last_idx*strides_[Dim];
    }

    template < typename... Ints >
    constexpr int index(Ints... idx) const {
        return index_part<0>(idx...);
    }

#ifdef ATLAS_ARRAYVIEW_BOUNDS_CHECKING
    template < typename... Ints >
    void check_bounds(Ints... idx) const {
        static_assert ( sizeof...(idx) == Rank , "Expected number of indices is different from rank of array" );
        return check_bounds_part<0>(idx...);
    }
#else
    template < typename... Ints >
    void check_bounds(Ints...) const {}
#endif

    template < typename... Ints >
    void check_bounds_force(Ints... idx) const {
        static_assert ( sizeof...(idx) == Rank , "Expected number of indices is different from rank of array" );
        return check_bounds_part<0>(idx...);
    }

    template < int Dim, typename Int, typename... Ints >
    void check_bounds_part(Int idx, Ints... next_idx) const {
        if( size_t(idx) >= shape_[Dim] ) {
            throw_OutOfRange( "LocalView", array_dim<Dim>(), idx, shape_[Dim] );
        }
        check_bounds_part<Dim+1>( next_idx... );
    }

    template < int Dim, typename Int >
    void check_bounds_part(Int last_idx) const {
        if( size_t(last_idx) >= shape_[Dim] ) {
            throw_OutOfRange( "LocalView", array_dim<Dim>(), last_idx, shape_[Dim] );
        }
    }

private:

// -- Private data

  value_type *data_;
  size_t size_;
  size_t shape_[Rank];
  size_t strides_[Rank];

};

} // namespace array
} // namespace atlas
