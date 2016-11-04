/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

/// @file ArrayView.h
/// This file contains the ArrayView class, a class that allows to wrap any contiguous raw data into
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
///     ArrayView<int,2> matrix( array, shape, strides );
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
/// There is also an easier way to wrap Field and Array classes:
///
/// Example 3:
///     ArrayView<int,3> fieldview( field::Field );
///     ArrayView<int,2> arrayview( Array );
///
/// @author Willem Deconinck
#ifndef atlas_ArrayView_h
#define atlas_ArrayView_h

#include <cstddef>
#include <vector>
#include "atlas/internals/atlas_defines.h"
#include "atlas/array/GridToolsTraits.h"

namespace atlas {
namespace array {
  class Array;
  template <typename DATA_TYPE, int RANK=0> class ArrayView;
} // namespace array
} // namespace atlas

//------------------------------------------------------------------------------------------------------

#include "atlas/array/ArrayUtil.h"
#include "atlas/array/ArrayView_iterator.h"

//------------------------------------------------------------------------------------------------------

namespace atlas {
namespace array {

#ifdef ATLAS_HAVE_GRIDTOOLS_STORAGE

template< typename DATA_TYPE, int RANK >
class ArrayView
{
public:
// -- Type definitions
  typedef ArrayView_iterator<DATA_TYPE,RANK>         iterator;
  typedef ArrayView_const_iterator<DATA_TYPE,RANK>   const_iterator;
  typedef typename remove_const<DATA_TYPE>::type  value_type;

  using data_view_t = data_view_tt<DATA_TYPE, RANK>;

public:

    ArrayView(data_view_t data_view) : gt_data_view_(data_view) {}

    template < typename... Coords, typename = typename boost::enable_if_c<(sizeof...(Coords) == RANK), int>::type >
    DATA_TYPE&
    GT_FUNCTION
    operator()(Coords... c) {
        assert(sizeof...(Coords) == RANK);
        return gt_data_view_(c...);
    }

    template <typename... Coords, typename = typename boost::enable_if_c<(sizeof...(Coords) == RANK), int>::type>
    GT_FUNCTION
    DATA_TYPE const& operator()(Coords... c) const {
      return gt_data_view_(c...);
    }

private:
    data_view_t gt_data_view_;
//// -- Constructors
//  ArrayView() {}
//  ArrayView( DATA_TYPE* data, const size_t size );
//  ArrayView( DATA_TYPE* data, const ArrayShape::value_type shape[1], const ArrayStrides::value_type strides[1] );
//  ArrayView( const DATA_TYPE* data, const ArrayShape::value_type shape[1] );
//  ArrayView( const DATA_TYPE* data, const ArrayShape& shape );
//  ArrayView( const Array& );

//// -- Iterators
//  iterator       begin();
//  iterator       end();
//  const_iterator begin() const;
//  const_iterator end()   const;
//  const_iterator cbegin() const;
//  const_iterator cend() const;

//// -- Operators
//  const DATA_TYPE& operator()(size_t i) const;
//        DATA_TYPE& operator()(size_t i);
//  const DATA_TYPE& operator()(const ArrayIdx& idx) const;
//        DATA_TYPE& operator()(const ArrayIdx& idx);
//  const DATA_TYPE& operator[](size_t i) const;
//        DATA_TYPE& operator[](size_t i);
//  const DATA_TYPE& at(size_t i) const;
//        DATA_TYPE& at(size_t i);
//  void operator=(const DATA_TYPE& scalar);

//// -- Accessors
//  const DATA_TYPE* data() const;
//        DATA_TYPE* data();
//  const ArrayStrides::value_type* strides() const;
//  const ArrayShape::value_type* shape() const;
//  ArrayShape::value_type shape(const size_t i) const;
//  ArrayStrides::value_type stride(size_t i) const;
//  size_t rank() const;
//  size_t size() const;

//private:
//// -- Private data
//  DATA_TYPE* data_;
//  ArrayStrides::value_type strides_[1];
//  ArrayShape::value_type   shape_[1];
};

#else


template< typename DATA_TYPE >
class ArrayView<DATA_TYPE,0>
{
public:

// -- Type definitions
  typedef ArrayView_iterator<DATA_TYPE,0>         iterator;
  typedef ArrayView_const_iterator<DATA_TYPE,0>   const_iterator;
  typedef typename remove_const<DATA_TYPE>::type  value_type;

public:

// -- Constructors
  ArrayView() {}
  ArrayView( const DATA_TYPE* data, const size_t rank, const ArrayShape::value_type shape[], const ArrayStrides::value_type strides[] );
  ArrayView( const Array& );

// -- Iterators
  iterator begin();
  iterator end();
  const_iterator begin() const;
  const_iterator end()   const;
  const_iterator cbegin() const;
  const_iterator cend() const;

// -- Operators
  const DATA_TYPE& operator()(size_t i) const;
        DATA_TYPE& operator()(size_t i);
  const DATA_TYPE& operator()(size_t i, size_t j) const;
        DATA_TYPE& operator()(size_t i, size_t j);
  const DATA_TYPE& operator()(size_t i, size_t j, size_t k) const;
        DATA_TYPE& operator()(size_t i, size_t j, size_t k);
  const DATA_TYPE& operator()(size_t i, size_t j, size_t k, size_t l) const;
        DATA_TYPE& operator()(size_t i, size_t j, size_t k, size_t l);
  const DATA_TYPE& operator()(size_t i, size_t j, size_t k, size_t l, size_t m) const;
        DATA_TYPE& operator()(size_t i, size_t j, size_t k, size_t l, size_t m);
  const DATA_TYPE& operator()(const ArrayIdx& idx) const;
        DATA_TYPE& operator()(const ArrayIdx& idx);
  void operator=(const DATA_TYPE& scalar);
  void resize(size_t size1, size_t size2);

// -- Accessors
  const DATA_TYPE* data() const;
        DATA_TYPE* data();
  const ArrayStrides::value_type* strides() const;
  const ArrayShape::value_type* shape() const;
  ArrayStrides::value_type stride(size_t i) const;
  ArrayShape::value_type shape(size_t i) const;
  size_t rank() const;
  size_t size() const;

private:
// -- Private data
  DATA_TYPE* data_;
  ArrayStrides strides_;
  ArrayShape shape_;  void resize(size_t size1, size_t size2);

  size_t rank_;
  size_t size_;
};

//------------------------------------------------------------------------------------------------------

template< typename DATA_TYPE >
class ArrayView < DATA_TYPE, 1 >
{
public:
// -- Type definitions
  typedef ArrayView_iterator<DATA_TYPE,1>         iterator;
  typedef ArrayView_const_iterator<DATA_TYPE,1>   const_iterator;
  typedef typename remove_const<DATA_TYPE>::type  value_type;

public:

// -- Constructors
  ArrayView() {}
  ArrayView( DATA_TYPE* data, const size_t size );
  ArrayView( DATA_TYPE* data, const ArrayShape::value_type shape[1], const ArrayStrides::value_type strides[1] );
  ArrayView( const DATA_TYPE* data, const ArrayShape::value_type shape[1] );
  ArrayView( const DATA_TYPE* data, const ArrayShape& shape );
  ArrayView( const Array& );

// -- Iterators
  iterator       begin();
  iterator       end();
  const_iterator begin() const;
  const_iterator end()   const;
  const_iterator cbegin() const;
  const_iterator cend() const;

// -- Operators
  const DATA_TYPE& operator()(size_t i) const;
        DATA_TYPE& operator()(size_t i);
  const DATA_TYPE& operator()(const ArrayIdx& idx) const;
        DATA_TYPE& operator()(const ArrayIdx& idx);
  const DATA_TYPE& operator[](size_t i) const;
        DATA_TYPE& operator[](size_t i);
  const DATA_TYPE& at(size_t i) const;
        DATA_TYPE& at(size_t i);
  void operator=(const DATA_TYPE& scalar);

// -- Accessors
  const DATA_TYPE* data() const;
        DATA_TYPE* data();
  const ArrayStrides::value_type* strides() const;
  const ArrayShape::value_type* shape() const;
  ArrayShape::value_type shape(const size_t i) const;
  ArrayStrides::value_type stride(size_t i) const;
  size_t rank() const;
  size_t size() const;

private:
// -- Private data
  DATA_TYPE* data_;
  ArrayStrides::value_type strides_[1];
  ArrayShape::value_type   shape_[1];
};

//------------------------------------------------------------------------------------------------------

template< typename DATA_TYPE >
class ArrayView < DATA_TYPE, 2 >
{  void resize(size_t size1, size_t size2);

public:
// -- Type definitions
  typedef ArrayView_iterator<DATA_TYPE,2>       iterator;
  typedef ArrayView_const_iterator<DATA_TYPE,2> const_iterator;
  typedef typename remove_const<DATA_TYPE>::type  value_type;

public:

// -- Constructors
  ArrayView() {}
  ArrayView( const DATA_TYPE* data, const ArrayShape::value_type shape[2], const ArrayStrides::value_type strides[2] );
  ArrayView( const DATA_TYPE* data, const ArrayShape::value_type shape[2] );
  ArrayView( const DATA_TYPE* data, const ArrayShape& shape );
  ArrayView( const Array& );

// -- Iterators
  iterator begin();
  iterator end();
  const_iterator begin() const;
  const_iterator end()   const;
  const_iterator cbegin() const;
  const_iterator cend() const;

// -- Operators
  const DATA_TYPE& operator()(size_t i, size_t j) const;
        DATA_TYPE& operator()(size_t i, size_t j);
  const DATA_TYPE& operator()(const ArrayIdx& idx) const;
        DATA_TYPE& operator()(const ArrayIdx& idx);
  const ArrayView<DATA_TYPE,1> operator[](size_t i) const;
        ArrayView<DATA_TYPE,1> operator[](size_t i);
  const ArrayView<DATA_TYPE,1> at(size_t i) const;
        ArrayView<DATA_TYPE,1> at(size_t i);
  void operator=(const DATA_TYPE& scalar);

// -- Accessors
  const DATA_TYPE* data() const;
        DATA_TYPE* data();
  const ArrayStrides::value_type* strides() const;
  const ArrayShape::value_type* shape() const;
  ArrayStrides::value_type stride(size_t i) const;
  ArrayShape::value_type shape(size_t i) const;
  size_t rank() const;
  size_t size() const;

private:
// -- Private data
  DATA_TYPE* data_;
  size_t size_;
  ArrayStrides::value_type strides_[2];
  ArrayShape::value_type shape_[2];
};

//------------------------------------------------------------------------------------------------------

template< typename DATA_TYPE >
class ArrayView < DATA_TYPE, 3 >
{
public:
// -- Type definitions
  typedef ArrayView_iterator<DATA_TYPE,3>         iterator;
  typedef ArrayView_const_iterator<DATA_TYPE,3>   const_iterator;
  typedef typename remove_const<DATA_TYPE>::type  value_type;

// -- Constructors
  ArrayView() {}
  ArrayView( const DATA_TYPE* data, const ArrayShape::value_type shape[3], const ArrayStrides::value_type strides [3] );
  ArrayView( const DATA_TYPE* data, const ArrayShape::value_type shape[3] );
  ArrayView( const DATA_TYPE* data, const ArrayShape& shape );
  ArrayView( const Array& );

// -- Iterators
  iterator begin();
  iterator end();
  const_iterator begin() const;
  const_iterator end()   const;
  const_iterator cbegin() const;
  const_iterator cend() const;

// -- Operators
  const DATA_TYPE& operator()(size_t i, size_t j, size_t k) const;
        DATA_TYPE& operator()(size_t i, size_t j, size_t k);
  const DATA_TYPE& operator()(const ArrayIdx& idx) const;
        DATA_TYPE& operator()(const ArrayIdx& idx);
  const ArrayView<DATA_TYPE,2> operator[](size_t i) const;
        ArrayView<DATA_TYPE,2> operator[](size_t i);
  const ArrayView<DATA_TYPE,2> at(size_t i) const;
        ArrayView<DATA_TYPE,2> at(size_t i);
  void operator=(const DATA_TYPE& scalar);

// -- Accessors
  const DATA_TYPE* data() const;
        DATA_TYPE* data();
  const ArrayStrides::value_type* strides() const;
  const ArrayShape::value_type* shape() const;
  ArrayStrides::value_type stride(size_t i) const;
  ArrayShape::value_type shape(size_t i) const;
  size_t rank() const;
  size_t size() const;

private:
// --Private data
  DATA_TYPE* data_;
  size_t size_;
  ArrayStrides::value_type strides_[3];
  ArrayShape::value_type shape_[3];
};

//------------------------------------------------------------------------------------------------------

template< typename DATA_TYPE >
class ArrayView < DATA_TYPE, 4 >
{
public:

// -- Type definitions
  typedef ArrayView_iterator<DATA_TYPE,4>       iterator;
  typedef ArrayView_const_iterator<DATA_TYPE,4> const_iterator;
  typedef typename remove_const<DATA_TYPE>::type  value_type;

// -- Constructors
  ArrayView() {}
  ArrayView( DATA_TYPE* data, const ArrayShape::value_type shape[4], const ArrayStrides::value_type strides[4] );
  ArrayView( const DATA_TYPE* data, const ArrayShape::value_type shape[4] );
  ArrayView( const DATA_TYPE* data, const ArrayShape& shape );
  ArrayView( const Array& );

// -- Iterators
  iterator begin();
  iterator end();
  const_iterator cbegin() const;
  const_iterator cend() const;
  const_iterator begin() const;
  const_iterator end()   const;

// -- Operators
  const DATA_TYPE& operator()(size_t i, size_t j, size_t k, size_t l) const;
        DATA_TYPE& operator()(size_t i, size_t j, size_t k, size_t l);
  const DATA_TYPE& operator()(const ArrayIdx& idx) const;
        DATA_TYPE& operator()(const ArrayIdx& idx);
  const ArrayView<DATA_TYPE,3> operator[](size_t i) const;
        ArrayView<DATA_TYPE,3> operator[](size_t i);
  void operator=(const DATA_TYPE& scalar);

// -- Accessors
  const DATA_TYPE* data() const;
        DATA_TYPE* data();
  const ArrayStrides::value_type* strides() const;
  const ArrayShape::value_type* shape() const;
  ArrayStrides::value_type stride(size_t i) const;
  ArrayShape::value_type shape(size_t i) const;
  size_t rank() const;
  size_t size() const;

private:
// -- Private data
  DATA_TYPE* data_;
  size_t size_;
  ArrayStrides::value_type strides_[4];
  ArrayShape::value_type shape_[4];
};

//------------------------------------------------------------------------------------------------------
#endif

} // namespace array
} // namespace atlas

#include "atlas/array/ArrayView_impl.h"

#endif
