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

#pragma once

#include "atlas/internals/atlas_defines.h"

#ifdef ATLAS_HAVE_GRIDTOOLS_STORAGE
#include "atlas/array/gridtools/GridToolsArrayView.h"
#else

#define atlas_ArrayView_h

#include <cstddef>
#include <cstring>
#include <cassert>
#include "atlas/array/ArrayUtil.h"
#include "atlas/array/LocalView.h"
#include "eckit/exception/Exceptions.h"

namespace atlas {
namespace array {
  class Array;
  template <typename DATA_TYPE, int RANK> class ArrayView;
} // namespace array
} // namespace atlas

//------------------------------------------------------------------------------------------------------

namespace atlas {
namespace array {

//------------------------------------------------------------------------------------------------------

template <typename Value, int Rank> class ArrayViewT {
public:

// -- Type definitions
  typedef typename remove_const<Value>::type  value_type;

  using degenerated_array_return_t = typename std::conditional<(Rank==1), value_type&, LocalView<value_type,Rank-1> >::type;

  template <typename ReturnType = degenerated_array_return_t, bool ToScalar = false>
  struct degenerate_local_array {
      degenerate_local_array(ArrayViewT<value_type, Rank> const& av) : av_(av) {}
      ArrayViewT<value_type, Rank> const& av_;
      ReturnType apply(const size_t i) const {
          return LocalView<value_type, Rank - 1>(av_.data_ + av_.strides_[0] * i, av_.shape_ + 1, av_.strides_ + 1);
      }
  };

  template <typename ReturnType>
  struct degenerate_local_array<ReturnType, true> {
      degenerate_local_array(ArrayViewT<value_type, Rank> const& av) : av_(av) {}
      ArrayViewT<value_type, Rank> const& av_;
      ReturnType apply(const size_t i) const {
          return *(av_.data_ + av_.strides_[0] * i);
      }
  };

public:

  // -- Constructors

    ArrayViewT( const value_type* data, const size_t shape[], const size_t strides[] ) :
        data_(const_cast<value_type*>(data)) {
        size_ = 1;
        for( size_t j=0; j<Rank; ++j ) {
            shape_[j] = shape[j];
            strides_[j] = strides[j];
            size_ *= shape_[j];
        }
    }

    ArrayViewT( const value_type* data, const size_t shape[] ) :
        data_(const_cast<value_type*>(data)) {
        size_ = 1;
        for( int j=Rank-1; j>=0; --j ) {
            shape_[j] = shape[j];
            strides_[j] = size_;
            size_ *= shape_[j];
        }
    }

    ArrayViewT( const value_type* data, const ArrayShape& shape ) :
        data_(const_cast<value_type*>(data)) {
        size_ = 1;
        for( int j=Rank-1; j>=0; --j ) {
            shape_[j]   = shape[j];
            strides_[j] = size_;
            size_ *= shape_[j];
        }
    }

// -- Access methods

    template < typename... Coords >
    value_type&
    operator()(Coords... c) {
        assert(sizeof...(Coords) == Rank);
        return data_[index(c...)];
    }

    template < typename... Coords >
    const value_type&
    operator()(Coords... c) const {
        assert(sizeof...(Coords) == Rank);
        return data_[index(c...)];
    }

    // TODO: implement
    //  const Value& operator()(const ArrayIdx& idx) const
    //        Value& operator()(const ArrayIdx& idx);

      const degenerated_array_return_t operator[](size_t i) const {
        return degenerate_local_array<degenerated_array_return_t, Rank==1>(*this).apply(i);
      }
      degenerated_array_return_t operator[](size_t i) {
        return degenerate_local_array<degenerated_array_return_t, Rank==1>(*this).apply(i);
      }

      const degenerated_array_return_t at(size_t i) const {
        if( i>= shape(0) ) throw eckit::OutOfRange(i,shape(0),Here());
        return degenerate_local_array<degenerated_array_return_t, Rank==1>(*this).apply(i);
      }

      degenerated_array_return_t at(size_t i) {
        if( i>= shape(0) ) throw eckit::OutOfRange(i,shape(0),Here());
        return degenerate_local_array<degenerated_array_return_t, Rank==1>(*this).apply(i);
      }

    size_t size() const { return size_;}
    size_t rank() const { return Rank; }
    const size_t* strides() const { return strides_; }
    const size_t* shape() const { return shape_; }
    size_t shape(size_t idx) const { return shape_[idx]; }
    size_t stride(size_t idx) const { return strides_[idx]; }
    value_type const* data() const { return data_; }
    value_type*       data()       { return data_; }
    bool contiguous() const { return (size_ == shape_[0]*strides_[0] ? true : false); }

    void assign(const value_type& value) {
        ASSERT( contiguous() );
        value_type* raw_data = data();
        for( size_t j=0; j<size_; ++j ) {
            raw_data[j] = value;
        }
    }

  void dump(std::ostream& os) const {
    ASSERT( contiguous() );
    const value_type* data_ = data();
    os << "size: " << size() << " , values: ";
    os << "[ ";
    for( size_t j=0; j<size(); ++ j )
      os << data_[j] << " ";
    os << "]";
  }

private:

// -- Private methods

    template < typename... Ints >
    constexpr int index_part(int cnt, int first, Ints... ints) const {
        return (cnt < Rank) ? first * strides_[cnt] + index_part(cnt + 1, ints..., first) : 0;
    }

    template < typename... Ints >
    constexpr int index(Ints... idx) const {
        return index_part(0, idx...);
    }

private:

// -- Private data
  value_type *data_;
  size_t size_;
  size_t shape_[Rank];
  size_t strides_[Rank];
};

#define ARRAYVIEW_NEWIMPL
#ifdef ARRAYVIEW_NEWIMPL

template< typename Value, int Rank >
class ArrayView : public ArrayViewT<Value,Rank> {
public:
  typedef typename remove_const<Value>::type  value_type;

  ArrayView( const value_type* data, const size_t shape[], const size_t strides[] ) :
    ArrayViewT<Value,Rank>(data,shape,strides) {}

  ArrayView( const value_type* data, const size_t shape[] ) :
    ArrayViewT<Value,Rank>(data,shape) {}

  ArrayView( const value_type* data, const ArrayShape& shape ) :
    ArrayViewT<Value,Rank>(data,shape) {}
};

} // namespace array
} // namespace atlas

#else

template< typename DATA_TYPE >
class ArrayView < DATA_TYPE, 1 >
{
public:
// -- Type definitions
//  typedef ArrayView_iterator<DATA_TYPE,1>         iterator;
//  typedef ArrayView_const_iterator<DATA_TYPE,1>   const_iterator;
  typedef typename remove_const<DATA_TYPE>::type  value_type;

public:

// -- Constructors
  ArrayView() {}
  ArrayView( DATA_TYPE* data, const size_t size );
  ArrayView( DATA_TYPE* data, const ArrayShape::value_type shape[1], const ArrayStrides::value_type strides[1] );
  ArrayView( const DATA_TYPE* data, const ArrayShape::value_type shape[1] );
  ArrayView( const DATA_TYPE* data, const ArrayShape& shape );
  // ArrayView( const Array& );

// -- Iterators
//  iterator       begin();
//  iterator       end();
//  const_iterator begin() const;
//  const_iterator end()   const;
//  const_iterator cbegin() const;
//  const_iterator cend() const;

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

  void assign(const DATA_TYPE& value) {
     ASSERT( contiguous() );
     DATA_TYPE* raw_data = data();
     for( size_t j=0; j<size(); ++j ) {
       raw_data[j] = value;
     }
  }

   bool contiguous() const
   {
     return true;
   }

  void dump(std::ostream& os) const
  {
    ASSERT( contiguous() );
    const DATA_TYPE* data_ = data();
    os << "size: " << size() << " , values: ";
    os << "[ ";
    for( size_t j=0; j<size(); ++ j )
      os << data_[j] << " ";
    os << "]";
  }


private:
// -- Private data
  DATA_TYPE* data_;
  ArrayStrides::value_type strides_[1];
  ArrayShape::value_type   shape_[1];
};

//------------------------------------------------------------------------------------------------------

template< typename DATA_TYPE >
class ArrayView < DATA_TYPE, 2 >
{
public:
// -- Type definitions
//  typedef ArrayView_iterator<DATA_TYPE,2>       iterator;
//  typedef ArrayView_const_iterator<DATA_TYPE,2> const_iterator;
  typedef typename remove_const<DATA_TYPE>::type  value_type;

public:

// -- Constructors
  ArrayView() {}
  ArrayView( const DATA_TYPE* data, const ArrayShape::value_type shape[2], const ArrayStrides::value_type strides[2] );
  ArrayView( const DATA_TYPE* data, const ArrayShape::value_type shape[2] );
  ArrayView( const DATA_TYPE* data, const ArrayShape& shape );
  // ArrayView( const Array& );

// -- Iterators
//  iterator begin();
//  iterator end();
//  const_iterator begin() const;
//  const_iterator end()   const;
//  const_iterator cbegin() const;
//  const_iterator cend() const;

// -- Operators
  const DATA_TYPE& operator()(size_t i, size_t j) const;
        DATA_TYPE& operator()(size_t i, size_t j);
  const DATA_TYPE& operator()(const ArrayIdx& idx) const;
        DATA_TYPE& operator()(const ArrayIdx& idx);
  const LocalView<DATA_TYPE,1> operator[](size_t i) const;
        LocalView<DATA_TYPE,1> operator[](size_t i);
  const LocalView<DATA_TYPE,1> at(size_t i) const;
        LocalView<DATA_TYPE,1> at(size_t i);
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

  void dump(std::ostream& os) const
  {
    ASSERT( contiguous() );
    const DATA_TYPE* data_ = data();
    os << "size: " << size() << " , values: ";
    os << "[ ";
    for( size_t j=0; j<size(); ++ j )
      os << data_[j] << " ";
    os << "]";
  }


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
//  typedef ArrayView_iterator<DATA_TYPE,3>         iterator;
//  typedef ArrayView_const_iterator<DATA_TYPE,3>   const_iterator;
  typedef typename remove_const<DATA_TYPE>::type  value_type;

// -- Constructors
  ArrayView() {}
  ArrayView( const DATA_TYPE* data, const ArrayShape::value_type shape[3], const ArrayStrides::value_type strides [3] );
  ArrayView( const DATA_TYPE* data, const ArrayShape::value_type shape[3] );
  ArrayView( const DATA_TYPE* data, const ArrayShape& shape );
  // ArrayView( const Array& );

// -- Iterators
//  iterator begin();
//  iterator end();
//  const_iterator begin() const;
//  const_iterator end()   const;
//  const_iterator cbegin() const;
//  const_iterator cend() const;

// -- Operators
  const DATA_TYPE& operator()(size_t i, size_t j, size_t k) const;
        DATA_TYPE& operator()(size_t i, size_t j, size_t k);
  const DATA_TYPE& operator()(const ArrayIdx& idx) const;
        DATA_TYPE& operator()(const ArrayIdx& idx);
  const LocalView<DATA_TYPE,2> operator[](size_t i) const;
        LocalView<DATA_TYPE,2> operator[](size_t i);
  const LocalView<DATA_TYPE,2> at(size_t i) const;
        LocalView<DATA_TYPE,2> at(size_t i);
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

  void dump(std::ostream& os) const
  {
    ASSERT( contiguous() );
    const DATA_TYPE* data_ = data();
    os << "size: " << size() << " , values: ";
    os << "[ ";
    for( size_t j=0; j<size(); ++ j )
      os << data_[j] << " ";
    os << "]";
  }

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
//  typedef ArrayView_iterator<DATA_TYPE,4>       iterator;
//  typedef ArrayView_const_iterator<DATA_TYPE,4> const_iterator;
  typedef typename remove_const<DATA_TYPE>::type  value_type;

// -- Constructors
  ArrayView() {}
  ArrayView( DATA_TYPE* data, const ArrayShape::value_type shape[4], const ArrayStrides::value_type strides[4] );
  ArrayView( const DATA_TYPE* data, const ArrayShape::value_type shape[4] );
  ArrayView( const DATA_TYPE* data, const ArrayShape& shape );
  // ArrayView( const Array& );

// -- Iterators
//  iterator begin();
//  iterator end();
//  const_iterator cbegin() const;
//  const_iterator cend() const;
//  const_iterator begin() const;
//  const_iterator end()   const;

// -- Operators
  const DATA_TYPE& operator()(size_t i, size_t j, size_t k, size_t l) const;
        DATA_TYPE& operator()(size_t i, size_t j, size_t k, size_t l);
  const DATA_TYPE& operator()(const ArrayIdx& idx) const;
        DATA_TYPE& operator()(const ArrayIdx& idx);
  const LocalView<DATA_TYPE,3> operator[](size_t i) const;
        LocalView<DATA_TYPE,3> operator[](size_t i);
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

  void dump(std::ostream& os) const
  {
    ASSERT( contiguous() );
    const DATA_TYPE* data_ = data();
    os << "size: " << size() << " , values: ";
    os << "[ ";
    for( size_t j=0; j<size(); ++ j )
      os << data_[j] << " ";
    os << "]";
  }

private:
// -- Private data
  DATA_TYPE* data_;
  size_t size_;
  ArrayStrides::value_type strides_[4];
  ArrayShape::value_type shape_[4];
};

//------------------------------------------------------------------------------------------------------

} // namespace array
} // namespace atlas

#include "atlas/array/native/ArrayView_impl.h"
#endif
#endif
