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
///     ArrayView<int,2> matrix( array, shape );
/// which is identical for this matrix to previous Example 1
///
/// @author Willem Deconinck
#pragma once

#include <cstddef>
#include <cstring>
#include <cassert>
#include <type_traits>
#include "atlas/internals/atlas_defines.h"
#include "atlas/array/ArrayUtil.h"

//------------------------------------------------------------------------------------------------------

namespace atlas {
namespace array {

  template< typename DATA_TYPE, int RANK >
  class LocalView
  {
  public:
    DATA_TYPE *data_;
    size_t shape_[RANK];
    size_t strides_[RANK];
    size_t size_;
  public:
    size_t size() const { return size_;}

  // -- Type definitions
    typedef typename remove_const<DATA_TYPE>::type  value_type;

    using degenerated_array_return_t = typename std::conditional<(RANK==1), DATA_TYPE&, LocalView<DATA_TYPE,RANK-1> >::type;

    template <typename ReturnType = degenerated_array_return_t, bool ToScalar = false>
    struct degenerate_local_array {
      degenerate_local_array(LocalView<DATA_TYPE, RANK> const& lv) : lv_(lv) {}
      LocalView<DATA_TYPE, RANK> const& lv_;
      ReturnType apply(const size_t i) const {
        return LocalView<DATA_TYPE, RANK - 1>(lv_.data_ + lv_.strides_[0] * i, lv_.shape_ + 1, lv_.strides_ + 1);
      }
    };

    template <typename ReturnType>
    struct degenerate_local_array<ReturnType, true> {
      degenerate_local_array(LocalView<DATA_TYPE, RANK> const& lv) : lv_(lv) {}

      LocalView<DATA_TYPE, RANK> const& lv_;
      ReturnType apply(const size_t i) const {
        return *(lv_.data_ + lv_.strides_[0] * i);
        ;
      }
    };

  public:

      LocalView( DATA_TYPE* data, const size_t shape[], const size_t strides[] ) :
        data_(data)
      {
        size_ = 1;
        for( size_t j=0; j<RANK; ++j ) {
          shape_[j] = shape[j];
          strides_[j] = strides[j];
          size_ *= shape_[j];
        }
      }

      LocalView( DATA_TYPE* data, const size_t shape[] ) :
        data_(data)
      {
        size_ = 1;
        for( int j=RANK-1; j>=0; --j ) {
          shape_[j] = shape[j];
          strides_[j] = size_;
          size_ *= shape_[j];
        }
      }

      LocalView( DATA_TYPE* data, const ArrayShape& shape ) :
        data_(data)
      {
        size_ = 1;
        for( int j=RANK-1; j>=0; --j ) {
          shape_[j]   = shape[j];
          strides_[j] = size_;
          size_ *= shape_[j];
        }
      }

private:
      template < typename... Ints >
      constexpr int index_part(int cnt, int first, Ints... ints) const {
          return (cnt < RANK) ? first * strides_[cnt] + index_part(cnt + 1, ints..., first) : 0;
      }

      template < typename... Ints >
      constexpr int index(Ints... idx) const {
          return index_part(0, idx...);
      }
public:
      template < typename... Coords >
      DATA_TYPE&
      operator()(Coords... c) {
          assert(sizeof...(Coords) == RANK);
          return data_[index(c...)];
      }

      template < typename... Coords >
      const DATA_TYPE&
      operator()(Coords... c) const {
          assert(sizeof...(Coords) == RANK);
          return data_[index(c...)];
      }

      size_t shape(size_t idx) const { return shape_[idx]; }

      degenerated_array_return_t at(const size_t i) const {
          return degenerate_local_array<degenerated_array_return_t, RANK==1>(*this).apply(i); }


      DATA_TYPE* data() { return data_; }
      DATA_TYPE const* data() const { return data_; }

      void assign(const DATA_TYPE& value) {
        ASSERT( contiguous() );
        DATA_TYPE* raw_data = data();
        for( size_t j=0; j<size_; ++j ) {
          raw_data[j] = value;
        }
      }

      bool contiguous() const
      {
        return (size_ == shape_[0]*strides_[0] ? true : false);
      }
  };

} // namespace array
} // namespace atlas
