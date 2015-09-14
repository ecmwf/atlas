/*
 * (C) Copyright 1996-2014 ECMWF.
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
/// Example:
/// int[] array = { 1, 2, 3, 4, 5, 6, 7, 8, 9};
/// int[2] strides = { 3, 1 };
/// int[2] shape = { 3, 3 };
/// ArrayView<int,2> matrix( array, strides, shape );
/// for( size_t i=0; i<matrix.shape(0); ++i ) {
///   for( size_t j=0; j<matrix.shape(1); ++j ) {
///
///  matrix(i,j) *= 10;
///   }
/// }
///
/// There is also an easier way to wrap Field and Array classes:
/// ArrayView<int,3> fieldview( Field );
/// ArrayView<int,2> arrayview( Array );

#ifndef atlas_ArrayView_h
#define atlas_ArrayView_h

#include <cstddef>
#include <vector>
#include <sstream>
#include <iterator>     // std::iterator, std::input_iterator_tag
#include <numeric> // std::accumulate
#include <functional> // std::multiplies
#define ATLAS_ARRAYVIEW_BOUNDS_CHECKING

#ifdef ATLAS_ARRAYVIEW_BOUNDS_CHECKING
#include <eckit/exception/Exceptions.h>


#define CHECK_RANK(R)\
  if(rank()!=R) { std::ostringstream msg; msg << "ArrayView  rank mismatch: rank()="<<rank()<< " != " << R; throw eckit::OutOfRange(msg.str(),Here()); }
#define CHECK_BOUNDS(idx) {\
  for( size_t d=0; d<rank(); ++d ) { \
    if(idx[d]>=shape_[d]) {std::ostringstream msg; msg << "index " << d << " out of bounds: " << idx[d] << " >= " << shape_[d]; throw eckit::OutOfRange(msg.str(),Here()); } } }
#define CHECK_BOUNDS_1(i)\
	if(i>=shape_[0]) {std::ostringstream msg; msg << "ArrayView(i) index out of bounds: i=" << i << " >= " << shape_[0]; throw eckit::OutOfRange(msg.str(),Here()); }
#define CHECK_BOUNDS_2(i,j)\
	if(i>=shape_[0]) {std::ostringstream msg; msg << "ArrayView(i,j) index out of bounds: i=" << i << " >= " << shape_[0]; throw eckit::OutOfRange(msg.str(),Here()); }\
	if(j>=shape_[1]) {std::ostringstream msg; msg << "ArrayView(i,j) index out of bounds: j=" << j << " >= " << shape_[1]; throw eckit::OutOfRange(msg.str(),Here()); }
#define CHECK_BOUNDS_3(i,j,k)\
	if(i>=shape_[0]) {std::ostringstream msg; msg << "ArrayView(i,j,k) index out of bounds: i=" << i << " >= " << shape_[0]; throw eckit::OutOfRange(msg.str(),Here()); }\
	if(j>=shape_[1]) {std::ostringstream msg; msg << "ArrayView(i,j,k) index out of bounds: j=" << j << " >= " << shape_[1]; throw eckit::OutOfRange(msg.str(),Here()); }\
	if(k>=shape_[2]) {std::ostringstream msg; msg << "ArrayView(i,j,k) index out of bounds: k=" << k << " >= " << shape_[2]; throw eckit::OutOfRange(msg.str(),Here()); }
#define CHECK_BOUNDS_4(i,j,k,l)\
	if(i>=shape_[0]) {std::ostringstream msg; msg << "ArrayView(i,j,k,l) index out of bounds: i=" << i << " >= " << shape_[0]; throw eckit::OutOfRange(msg.str(),Here()); }\
	if(j>=shape_[1]) {std::ostringstream msg; msg << "ArrayView(i,j,k,l) index out of bounds: j=" << j << " >= " << shape_[1]; throw eckit::OutOfRange(msg.str(),Here()); }\
	if(k>=shape_[2]) {std::ostringstream msg; msg << "ArrayView(i,j,k,l) index out of bounds: k=" << k << " >= " << shape_[2]; throw eckit::OutOfRange(msg.str(),Here()); }\
	if(l>=shape_[3]) {std::ostringstream msg; msg << "ArrayView(i,j,k,l) index out of bounds: l=" << l << " >= " << shape_[3]; throw eckit::OutOfRange(msg.str(),Here()); }
#define CHECK_BOUNDS_5(i,j,k,l,m)\
	if(i>=shape_[0]) {std::ostringstream msg; msg << "ArrayView(i,j,k,l,m) index out of bounds: i=" << i << " >= " << shape_[0]; throw eckit::OutOfRange(msg.str(),Here()); }\
	if(j>=shape_[1]) {std::ostringstream msg; msg << "ArrayView(i,j,k,l,m) index out of bounds: j=" << j << " >= " << shape_[1]; throw eckit::OutOfRange(msg.str(),Here()); }\
	if(k>=shape_[2]) {std::ostringstream msg; msg << "ArrayView(i,j,k,l,m) index out of bounds: k=" << k << " >= " << shape_[2]; throw eckit::OutOfRange(msg.str(),Here()); }\
	if(l>=shape_[3]) {std::ostringstream msg; msg << "ArrayView(i,j,k,l,m) index out of bounds: l=" << l << " >= " << shape_[3]; throw eckit::OutOfRange(msg.str(),Here()); }\
	if(m>=shape_[4]) {std::ostringstream msg; msg << "ArrayView(i,j,k,l,m) index out of bounds: m=" << m << " >= " << shape_[4]; throw eckit::OutOfRange(msg.str(),Here()); }
#else
#define CHECK_RANK(R)
#define CHECK_BOUNDS(i,max)
#define CHECK_BOUNDS_1(i)
#define CHECK_BOUNDS_2(i,j)
#define CHECK_BOUNDS_3(i,j,k)
#define CHECK_BOUNDS_4(i,j,k,l)
#define CHECK_BOUNDS_5(i,j,k,l,m)
#endif

#include "ArrayUtil.h"


//------------------------------------------------------------------------------------------------------

namespace atlas {

class Array;
class Field;
template< typename DATA_TYPE > class ArrayT;

template <typename DATA_TYPE, int RANK=0> class ArrayView;
template <typename DATA_TYPE, int RANK=0> class ArrayView_iterator;
template <typename DATA_TYPE, int RANK=0> class ArrayView_const_iterator;

template <typename DATA_TYPE, int RANK>
class ArrayView_iterator : public std::iterator<std::forward_iterator_tag, DATA_TYPE>
{
    friend class ArrayView_const_iterator<DATA_TYPE,RANK>;
private:
  DATA_TYPE* p_;                    //!< raw pointer
  ArrayView<DATA_TYPE,RANK>* arr_;  //!< array to be iterated
  ArrayIdx loc_;                    //!< current position in array
  int fastest_idx_;                 //!< store fastest-moving index

public:

  /// @brief constructor (to be used for end of range)
  ArrayView_iterator();

  /// @brief constructor (to be used for begin of range)
  ArrayView_iterator(ArrayView<DATA_TYPE,RANK>* arr);

  /// @brief constructor from array at specific point
  ArrayView_iterator(ArrayView<DATA_TYPE,RANK>* arr, const ArrayIdx& loc);

  /// @brief constructor from other iterator
  ArrayView_iterator(const ArrayView_iterator& it);

  /// @brief pre-increment operator
  ArrayView_iterator& operator++()    { return increment(fastest_idx_); }

  /// @brief post-increment operator
  ArrayView_iterator operator++(int)  { ArrayView_iterator<DATA_TYPE,RANK> tmp(*this); operator++(); return tmp; }

  /// @brief equals operator
  bool operator==(const ArrayView_iterator& rhs) { return p_==rhs.p_; }

  /// @brief not-equals operator
  bool operator!=(const ArrayView_iterator& rhs) { return p_!=rhs.p_; }

  /// @brief dereference operator
  DATA_TYPE& operator*()      { return *p_; }

  /// @brief current position in array
  const ArrayIdx& pos() const { return loc_; }

private:
  ArrayView_iterator& increment(int d);

};

template <typename DATA_TYPE, int RANK>
class ArrayView_const_iterator : public std::iterator<std::forward_iterator_tag, DATA_TYPE>
{
  friend class ArrayView_iterator<DATA_TYPE,RANK>;
private:
  const DATA_TYPE* p_;                    //!< raw pointer
  const ArrayView<DATA_TYPE,RANK>* arr_;  //!< array to be iterated
  ArrayIdx loc_;                    //!< current position in array
  int fastest_idx_;                 //!< store fastest-moving index

public:

  /// @brief constructor (to be used for end of range)
  ArrayView_const_iterator();

  /// @brief constructor (to be used for begin of range)
  ArrayView_const_iterator(const ArrayView<DATA_TYPE,RANK>* arr);

  /// @brief constructor from array at specific point
  ArrayView_const_iterator(const ArrayView<DATA_TYPE,RANK>* arr, const ArrayIdx& loc);

  /// @brief constructor from other iterator
  ArrayView_const_iterator(const ArrayView_iterator<DATA_TYPE,RANK>& it);

  /// @brief constructor from other iterator
  ArrayView_const_iterator(const ArrayView_const_iterator& it);

  /// @brief pre-increment operator
  ArrayView_const_iterator& operator++()    { return increment(fastest_idx_); }

  /// @brief post-increment operator
  ArrayView_const_iterator operator++(int)  { ArrayView_const_iterator<DATA_TYPE,RANK> tmp(*this); operator++(); return tmp; }

  /// @brief equals operator
  bool operator==(const ArrayView_iterator<DATA_TYPE,RANK>& rhs) { return p_==rhs.p_; }

  /// @brief equals operator
  bool operator==(const ArrayView_const_iterator& rhs) { return p_==rhs.p_; }

  /// @brief not-equals operator
  bool operator!=(const ArrayView_iterator<DATA_TYPE,RANK>& rhs) { return p_!=rhs.p_; }

  /// @brief not-equals operator
  bool operator!=(const ArrayView_const_iterator& rhs) { return p_!=rhs.p_; }

  /// @brief dereference operator
  const DATA_TYPE& operator*() const { return *p_; }

  /// @brief current position in array
  const ArrayIdx& pos() const { return loc_; }

private:
  ArrayView_const_iterator& increment(int d);
};

//------------------------------------------------------------------------------------------------------

template< typename DATA_TYPE >
class ArrayView<DATA_TYPE,0>
{
public:
  typedef ArrayView_iterator<DATA_TYPE,0>       iterator;
  typedef ArrayView_const_iterator<DATA_TYPE,0> const_iterator;
  typedef typename remove_const<DATA_TYPE>::type  value_type;
public:
  ArrayView() {}
  ArrayView( const DATA_TYPE* data,
             const ArrayStrides::value_type strides[],
             const ArrayShape::value_type shape[],
             const size_t rank ):
    data_( const_cast<DATA_TYPE*>(data) ), rank_(rank)
  {
    strides_.assign(strides,strides+rank_);
    shape_.assign(shape,shape+rank_);
    size_ = std::accumulate(shape_.data(),shape_.data()+rank_,1,std::multiplies<size_t>());
  }

  ArrayView( const ArrayT<DATA_TYPE>& );
  ArrayView( const Array& );
  ArrayView( const Field& );
  iterator begin() { return iterator(this); }
  iterator end()   { return iterator(); }
  const_iterator begin() const  { return const_iterator(this); }
  const_iterator end()   const  { return const_iterator(); }
  const DATA_TYPE& operator()(size_t i) const { CHECK_RANK(1); CHECK_BOUNDS_1(i); return *(data_+strides_[0]*i); }
  DATA_TYPE&       operator()(size_t i)       { CHECK_RANK(1); CHECK_BOUNDS_1(i); return *(data_+strides_[0]*i); }
  const DATA_TYPE& operator()(size_t i, size_t j) const  { CHECK_RANK(2); CHECK_BOUNDS_2(i,j); return *(data_+strides_[0]*i+j*strides_[1]); }
  DATA_TYPE&       operator()(size_t i, size_t j)        { CHECK_RANK(2); CHECK_BOUNDS_2(i,j); return *(data_+strides_[0]*i+j*strides_[1]); }
  const DATA_TYPE& operator()(size_t i, size_t j, size_t k) const { CHECK_RANK(3); CHECK_BOUNDS_3(i,j,k); return *(data_+strides_[0]*i+j*strides_[1]+k*strides_[2]); }
  DATA_TYPE&       operator()(size_t i, size_t j, size_t k)       { CHECK_RANK(3); CHECK_BOUNDS_3(i,j,k); return *(data_+strides_[0]*i+j*strides_[1]+k*strides_[2]); }
  const DATA_TYPE& operator()(size_t i, size_t j, size_t k, size_t l) const { CHECK_RANK(4); CHECK_BOUNDS_4(i,j,k,l); return *(data_+strides_[0]*i+j*strides_[1]+k*strides_[2]+l*strides_[3]); }
  DATA_TYPE&       operator()(size_t i, size_t j, size_t k, size_t l)       { CHECK_RANK(4); CHECK_BOUNDS_4(i,j,k,l); return *(data_+strides_[0]*i+j*strides_[1]+k*strides_[2]+l*strides_[3]); }
  const DATA_TYPE& operator()(size_t i, size_t j, size_t k, size_t l, size_t m) const { CHECK_RANK(5); CHECK_BOUNDS_5(i,j,k,l,m); return *(data_+strides_[0]*i+j*strides_[1]+k*strides_[2]+l*strides_[3]+m*strides_[4]); }
  DATA_TYPE&       operator()(size_t i, size_t j, size_t k, size_t l, size_t m)       { CHECK_RANK(5); CHECK_BOUNDS_5(i,j,k,l,m); return *(data_+strides_[0]*i+j*strides_[1]+k*strides_[2]+l*strides_[3]+m*strides_[4]); }
  DATA_TYPE& operator()(const ArrayIdx& idx)
  {
    CHECK_RANK(idx.size()); CHECK_BOUNDS(idx);
    size_t p=0;
    for( size_t d=0; d<rank(); ++d )
      p += idx[d]*strides_[d];
    return *(data_+p);
  }
  const DATA_TYPE& operator()(const ArrayIdx& idx) const
  {
    CHECK_RANK(idx.size()); CHECK_BOUNDS(idx);
    size_t p=0;
    for( size_t d=0; d<rank(); ++d )
      p += idx[d]*strides_[d];
    return *(data_+p);
  }
  DATA_TYPE*       data()        { return data_; }
  const DATA_TYPE* data() const  { return data_; }
  const ArrayStrides::value_type* strides() const  { return strides_.data(); }
  const ArrayShape::value_type* shape() const  { return shape_.data(); }
  ArrayStrides::value_type stride(size_t i) const { return strides_[i]; }
  ArrayShape::value_type shape(size_t i) const { return shape_[i]; }
  size_t rank() const       { return rank_; }
  size_t size() const       { return size_; }
  void operator=(const DATA_TYPE& scalar) { for(size_t n=0; n<size_; ++n) *(data_+n)=scalar; }
private:
  DATA_TYPE* data_;
  ArrayStrides strides_;
  ArrayShape shape_;
  size_t rank_;
  size_t size_;
};


//------------------------------------------------------------------------------------------------------

template< typename DATA_TYPE >
class ArrayView < DATA_TYPE, 1 >
{
public:
  typedef ArrayView_iterator<DATA_TYPE,1>       iterator;
  typedef ArrayView_const_iterator<DATA_TYPE,1> const_iterator;
  typedef typename remove_const<DATA_TYPE>::type  value_type;
public:
  ArrayView() {}
  ArrayView( DATA_TYPE* data, const size_t size ) : data_( const_cast<DATA_TYPE*>(data) )
  {
    shape_[0]=size; strides_[0]=1;
  }
  ArrayView( DATA_TYPE* data, const ArrayStrides::value_type strides[1], const ArrayShape::value_type shape[1] ) : data_( const_cast<DATA_TYPE*>(data) )
  {
    strides_[0]=strides[0];       shape_[0]=shape[0];
  }
  ArrayView( const DATA_TYPE* data, const ArrayShape::value_type shape[1] ) : data_( const_cast<DATA_TYPE*>(data) )
  {
    shape_[0]=shape[0]; strides_[0]=1;
  }
  ArrayView( const DATA_TYPE* data, const std::vector<size_t>& shape ) : data_( const_cast<DATA_TYPE*>(data) )
  {
    shape_[0]=shape[0]; strides_[0]=1;
  }
  ArrayView( const ArrayT<DATA_TYPE>& );
  ArrayView( const Array& );
  ArrayView( const Field& );

  iterator begin() { return iterator(this); }
  iterator end()   { return iterator(); }
  const_iterator begin() const  { return const_iterator(this); }
  const_iterator end()   const  { return const_iterator(); }
  const DATA_TYPE& operator()(size_t i) const { CHECK_BOUNDS_1(i); return *(data_+strides_[0]*i); }
  DATA_TYPE&       operator()(size_t i)       { CHECK_BOUNDS_1(i); return *(data_+strides_[0]*i); }
  const DATA_TYPE& operator()(const ArrayIdx& idx) const
  {
    CHECK_RANK(idx.size()); CHECK_BOUNDS(idx);
    size_t p=0;
    for( size_t d=0; d<rank(); ++d )
      p += idx[d]*strides_[d];
    return *(data_+p);
  }
  DATA_TYPE& operator()(const ArrayIdx& idx)
  {
    CHECK_RANK(idx.size()); CHECK_BOUNDS(idx);
    size_t p=0;
    for( size_t d=0; d<rank(); ++d )
      p += idx[d]*strides_[d];
    return *(data_+p);
  }
  const DATA_TYPE& operator[](size_t i) const { CHECK_BOUNDS_1(i); return *(data_+strides_[0]*i); }
  DATA_TYPE&       operator[](size_t i)       { CHECK_BOUNDS_1(i); return *(data_+strides_[0]*i); }
  DATA_TYPE*       data()        { return data_; }
  const DATA_TYPE* data() const  { return data_; }
  const ArrayStrides::value_type* strides() const   { return strides_; }
  const ArrayShape::value_type* shape() const   { return shape_; }
	ArrayShape::value_type shape(const size_t i) const { return shape_[0]; }
  ArrayStrides::value_type stride(size_t i) const { return strides_[0]; }
  size_t rank() const { return 1; }
  size_t size() const { return shape_[0]; }
  size_t total_size() const { return shape_[0]; }
  void operator=(const DATA_TYPE& scalar) { for(size_t n=0; n<total_size(); ++n) *(data_+n)=scalar; }
private:
  DATA_TYPE* data_;
  ArrayStrides::value_type strides_[1];
  ArrayShape::value_type   shape_[1];
};

//------------------------------------------------------------------------------------------------------

template< typename DATA_TYPE >
class ArrayView < DATA_TYPE, 2 >
{
public:
  typedef ArrayView_iterator<DATA_TYPE,2>       iterator;
  typedef ArrayView_const_iterator<DATA_TYPE,2> const_iterator;
  typedef typename remove_const<DATA_TYPE>::type  value_type;
public:

  ArrayView() {}
  ArrayView( const DATA_TYPE* data, const ArrayStrides::value_type strides[2], const ArrayShape::value_type shape[2] ) : data_( const_cast<DATA_TYPE*>(data) )
  {
    strides_[0]=strides[0];            shape_[0]=shape[0];
    strides_[1]=strides[1];            shape_[1]=shape[1];
  }
  ArrayView( const DATA_TYPE* data, const ArrayShape::value_type shape[2] ) : data_( const_cast<DATA_TYPE*>(data) )
  {
    shape_[0]=shape[0]; strides_[0]=shape[1];
    shape_[1]=shape[1]; strides_[1]=1;
  }
  ArrayView( const DATA_TYPE* data, const std::vector<size_t>& shape ) : data_( const_cast<DATA_TYPE*>(data) )
  {
    shape_[0]=shape[0]; strides_[0]=shape[1];
    shape_[1]=shape[1]; strides_[1]=1;
  }
  ArrayView( const ArrayT<DATA_TYPE>& );
  ArrayView( const Array& );
  ArrayView( const Field& );

  iterator begin() { return iterator(this); }
  iterator end()   { return iterator(); }
  const_iterator begin() const  { return const_iterator(this); }
  const_iterator end()   const  { return const_iterator(); }
  const DATA_TYPE& operator()(size_t i, size_t j) const  { CHECK_BOUNDS_2(i,j); return *(data_+strides_[0]*i+j*strides_[1]); }
  DATA_TYPE&       operator()(size_t i, size_t j)        { CHECK_BOUNDS_2(i,j); return *(data_+strides_[0]*i+j*strides_[1]); }
  const DATA_TYPE& operator()(const ArrayIdx& idx) const
  {
    CHECK_RANK(idx.size()); CHECK_BOUNDS(idx);
    size_t p=0;
    for( size_t d=0; d<rank(); ++d )
      p += idx[d]*strides_[d];
    return *(data_+p);
  }
  DATA_TYPE& operator()(const ArrayIdx& idx)
  {
    CHECK_RANK(idx.size()); CHECK_BOUNDS(idx);
    size_t p=0;
    for( size_t d=0; d<rank(); ++d )
      p += idx[d]*strides_[d];
    return *(data_+p);
  }
  const ArrayView<DATA_TYPE,1> operator[](size_t i) const { CHECK_BOUNDS_1(i); return ArrayView<DATA_TYPE,1>( data_+strides_[0]*i, strides_+1, shape_+1 ); }
  ArrayView<DATA_TYPE,1>       operator[](size_t i)       { CHECK_BOUNDS_1(i); return ArrayView<DATA_TYPE,1>( data_+strides_[0]*i, strides_+1, shape_+1 ); }
  DATA_TYPE*       data()            { return data_; }
  const DATA_TYPE* data() const      { return data_; }
  const ArrayStrides::value_type* strides() const   { return strides_; }
  const ArrayShape::value_type* shape() const   { return shape_; }
  ArrayStrides::value_type stride(size_t i) const { return strides_[i]; }
  ArrayShape::value_type shape(size_t i) const {return shape_[i]; }
  size_t rank() const { return 2; }
  size_t size() const { return shape_[0]; }
  size_t total_size() const { return shape_[0]*shape_[1]; }
  void operator=(const DATA_TYPE& scalar) { for(size_t n=0; n<total_size(); ++n) *(data_+n)=scalar; }
private:
  DATA_TYPE* data_;
  ArrayStrides::value_type strides_[2];
  ArrayShape::value_type shape_[2];
};

//------------------------------------------------------------------------------------------------------

template< typename DATA_TYPE >
class ArrayView < DATA_TYPE, 3 >
{
public:
  typedef ArrayView_iterator<DATA_TYPE,3>       iterator;
  typedef ArrayView_const_iterator<DATA_TYPE,3> const_iterator;
  typedef typename remove_const<DATA_TYPE>::type  value_type;
public:
  ArrayView() {}
  ArrayView( const DATA_TYPE* data, const ArrayStrides::value_type strides [3], const ArrayShape::value_type shape[3] ) : data_( const_cast<DATA_TYPE*>(data) )
  {
    strides_[0]=strides[0];            shape_[0]=shape[0];
    strides_[1]=strides[1];            shape_[1]=shape[1];
    strides_[2]=strides[2];            shape_[2]=shape[2];
  }
  ArrayView( const DATA_TYPE* data, const ArrayShape::value_type shape[3] ) : data_( const_cast<DATA_TYPE*>(data) )
  {
    shape_[0]=shape[0]; strides_[0]=shape[2]*shape[1];
    shape_[1]=shape[1]; strides_[1]=shape[2];
    shape_[2]=shape[2]; strides_[2]=1;
  }
  ArrayView( const DATA_TYPE* data, const std::vector<size_t>& shape ) : data_( const_cast<DATA_TYPE*>(data) )
  {
    shape_[0]=shape[0]; strides_[0]=shape[2]*shape[1];
    shape_[1]=shape[1]; strides_[1]=shape[2];
    shape_[2]=shape[2]; strides_[2]=1;
  }
  ArrayView( const ArrayT<DATA_TYPE>& );
  ArrayView( const Array& );
  ArrayView( const Field& );

  iterator begin() { return iterator(this); }
  iterator end()   { return iterator(); }
  const_iterator begin() const  { return const_iterator(this); }
  const_iterator end()   const  { return const_iterator(); }
  const DATA_TYPE& operator()(size_t i, size_t j, size_t k) const { CHECK_BOUNDS_3(i,j,k); return *(data_+strides_[0]*i+j*strides_[1]+k*strides_[2]); }
  DATA_TYPE&       operator()(size_t i, size_t j, size_t k)       { CHECK_BOUNDS_3(i,j,k); return *(data_+strides_[0]*i+j*strides_[1]+k*strides_[2]); }
  const DATA_TYPE& operator()(const ArrayIdx& idx) const
  {
    CHECK_RANK(idx.size()); CHECK_BOUNDS(idx);
    size_t p=0;
    for( size_t d=0; d<rank(); ++d )
      p += idx[d]*strides_[d];
    return *(data_+p);
  }
  DATA_TYPE& operator()(const ArrayIdx& idx)
  {
    CHECK_RANK(idx.size()); CHECK_BOUNDS(idx);
    size_t p=0;
    for( size_t d=0; d<rank(); ++d )
      p += idx[d]*strides_[d];
    return *(data_+p);
  }
  const ArrayView<DATA_TYPE,2> operator[](size_t i) const { CHECK_BOUNDS_1(i); return ArrayView<DATA_TYPE,2>( data_+strides_[0]*i, strides_+1, shape_+1 ); }
  ArrayView<DATA_TYPE,2>       operator[](size_t i)       { CHECK_BOUNDS_1(i);return ArrayView<DATA_TYPE,2>( data_+strides_[0]*i, strides_+1, shape_+1 ); }
  DATA_TYPE*       data()            { return data_; }
  const DATA_TYPE* data() const      { return data_; }
  const ArrayStrides::value_type* strides() const   { return strides_; }
  const ArrayShape::value_type* shape() const   { return shape_; }
  ArrayStrides::value_type stride(size_t i) const { return strides_[i]; }
  ArrayShape::value_type shape(size_t i) const { return shape_[i]; }
  size_t rank() const { return 3; }
  size_t size() const { return shape_[0]; }
  size_t total_size() const { return shape_[0]*shape_[1]*shape_[2]; }
  void operator=(const DATA_TYPE& scalar) { for(size_t n=0; n<total_size(); ++n) *(data_+n)=scalar; }
private:
  DATA_TYPE* data_;
  ArrayStrides::value_type strides_[3];
  ArrayShape::value_type shape_[3];
};

//------------------------------------------------------------------------------------------------------

template< typename DATA_TYPE >
class ArrayView < DATA_TYPE, 4 >
{
public:
  typedef ArrayView_iterator<DATA_TYPE,4>       iterator;
  typedef ArrayView_const_iterator<DATA_TYPE,4> const_iterator;
  typedef typename remove_const<DATA_TYPE>::type  value_type;
public:
  ArrayView() {}
  ArrayView( DATA_TYPE* data, const ArrayStrides::value_type strides[4], const ArrayShape::value_type shape[4] ) : data_( const_cast<DATA_TYPE*>(data) )
  {
    strides_[0]=strides[0];            shape_[0]=shape[0];
    strides_[1]=strides[1];            shape_[1]=shape[1];
    strides_[2]=strides[2];            shape_[2]=shape[2];
    strides_[3]=strides[3];            shape_[3]=shape[3];
  }
  ArrayView( const DATA_TYPE* data, const std::vector<size_t>& shape ) : data_( const_cast<DATA_TYPE*>(data) )
  {
    shape_[0]=shape[0]; strides_[0]=shape[3]*shape[2]*shape[1];
    shape_[1]=shape[1]; strides_[1]=shape[3]*shape[2];
    shape_[2]=shape[2]; strides_[2]=shape[3];
    shape_[3]=shape[3]; strides_[3]=1;
  }
  ArrayView( const DATA_TYPE* data, const ArrayShape::value_type shape[4] ) : data_( const_cast<DATA_TYPE*>(data) )
  {
    shape_[0]=shape[0]; strides_[0]=shape[3]*shape[2]*shape[1];
    shape_[1]=shape[1]; strides_[1]=shape[3]*shape[2];
    shape_[2]=shape[2]; strides_[2]=shape[3];
    shape_[3]=shape[3]; strides_[3]=1;
  }
  ArrayView( const ArrayT<DATA_TYPE>& );
  ArrayView( const Array& );
  ArrayView( const Field& );

  iterator begin() { return iterator(this); }
  iterator end()   { return iterator(); }
  const_iterator begin() const  { return const_iterator(this); }
  const_iterator end()   const  { return const_iterator(); }
  const DATA_TYPE& operator()(size_t i, size_t j, size_t k, size_t l) const { CHECK_BOUNDS_4(i,j,k,l); return *(data_+strides_[0]*i+j*strides_[1]+k*strides_[2]+l*strides_[3]); }
  DATA_TYPE&       operator()(size_t i, size_t j, size_t k, size_t l)       { CHECK_BOUNDS_4(i,j,k,l); return *(data_+strides_[0]*i+j*strides_[1]+k*strides_[2]+l*strides_[3]); }
  const DATA_TYPE& operator()(const ArrayIdx& idx) const
  {
    CHECK_RANK(idx.size()); CHECK_BOUNDS(idx);
    size_t p=0;
    for( size_t d=0; d<rank(); ++d )
      p += idx[d]*strides_[d];
    return *(data_+p);
  }
  DATA_TYPE& operator()(const ArrayIdx& idx)
  {
    CHECK_RANK(idx.size()); CHECK_BOUNDS(idx);
    size_t p=0;
    for( size_t d=0; d<rank(); ++d )
      p += idx[d]*strides_[d];
    return *(data_+p);
  }
  const ArrayView<DATA_TYPE,3> operator[](size_t i) const { CHECK_BOUNDS_1(i); return ArrayView<DATA_TYPE,3>( data_+strides_[0]*i, strides_+1, shape_+1 ); }
  ArrayView<DATA_TYPE,3>       operator[](size_t i)       { CHECK_BOUNDS_1(i); return ArrayView<DATA_TYPE,3>( data_+strides_[0]*i, strides_+1, shape_+1 ); }
  DATA_TYPE*       data()            { return data_; }
  const DATA_TYPE* data() const      { return data_; }
  const ArrayStrides::value_type* strides() const   { return strides_; }
  const ArrayShape::value_type* shape() const   { return shape_; }
  ArrayStrides::value_type stride(size_t i) const { return strides_[i]; }
  ArrayShape::value_type shape(size_t i) const { return shape_[i]; }
  size_t rank() const { return 4; }
  size_t size() const { return shape_[0]; }
  size_t total_size() const { return shape_[0]*shape_[1]*shape_[2]*shape_[3]; }
  void operator=(const DATA_TYPE& scalar) { for(size_t n=0; n<total_size(); ++n) *(data_+n)=scalar; }

private:
  DATA_TYPE* data_;
  ArrayStrides::value_type strides_[4];
  ArrayShape::value_type shape_[4];
};

//------------------------------------------------------------------------------------------------------

template <typename DATA_TYPE, int RANK>
ArrayView_iterator<DATA_TYPE,RANK>::ArrayView_iterator()
    : p_(0), arr_(0), loc_(0), fastest_idx_(0)
{
}

template <typename DATA_TYPE, int RANK>
ArrayView_iterator<DATA_TYPE,RANK>::ArrayView_iterator(ArrayView<DATA_TYPE,RANK>* arr):
 arr_(arr),loc_(arr->rank(),0)
{
  fastest_idx_ = arr_->rank()-1;
  p_ = arr_->data();
}

template <typename DATA_TYPE, int RANK>
ArrayView_iterator<DATA_TYPE,RANK>::ArrayView_iterator(ArrayView<DATA_TYPE,RANK>* arr, const ArrayIdx& loc):
 arr_(arr),loc_(loc)
{
  fastest_idx_ = arr_->rank()-1;
  p_ = &arr_->operator()(loc_);
}

template <typename DATA_TYPE, int RANK>
ArrayView_iterator<DATA_TYPE,RANK>::ArrayView_iterator(const ArrayView_iterator& it)
    : p_(it.p_), arr_(it.arr_), loc_(it.loc_)
{}

template <typename DATA_TYPE, int RANK>
ArrayView_iterator<DATA_TYPE,RANK>& ArrayView_iterator<DATA_TYPE,RANK>::increment(int d)
{
  ++loc_[d];
  p_ += arr_->stride(d);
  if( loc_[d] == arr_->shape(d) )
  {
    p_ -= arr_->stride(d)*arr_->shape(d);
    loc_[d]=0;
    if(d==0)
    {
      p_ = 0;
      return *this;
    }
    else
    {
      return increment(d-1);
    }
  }
  return *this;
}

//------------------------------------------------------------------------------------------------------

template <typename DATA_TYPE, int RANK>
ArrayView_const_iterator<DATA_TYPE,RANK>::ArrayView_const_iterator()
    : p_(0), arr_(0), loc_(0), fastest_idx_(0)
{
}

template <typename DATA_TYPE, int RANK>
ArrayView_const_iterator<DATA_TYPE,RANK>::ArrayView_const_iterator(const ArrayView<DATA_TYPE,RANK>* arr):
 arr_(arr),loc_(arr->rank(),0)
{
  fastest_idx_ = arr_->rank()-1;
  p_ = arr_->data();
}


template <typename DATA_TYPE, int RANK>
ArrayView_const_iterator<DATA_TYPE,RANK>::ArrayView_const_iterator(const ArrayView<DATA_TYPE,RANK>* arr, const ArrayIdx& loc):
 arr_(arr),loc_(loc)
{
  fastest_idx_ = arr_->rank()-1;
  p_ = &arr_->operator()(loc_);
}

template <typename DATA_TYPE, int RANK>
ArrayView_const_iterator<DATA_TYPE,RANK>::ArrayView_const_iterator(const ArrayView_iterator<DATA_TYPE,RANK>& it)
    : p_(it.p_), arr_(it.arr_), loc_(it.loc_)
{}

template <typename DATA_TYPE, int RANK>
ArrayView_const_iterator<DATA_TYPE,RANK>::ArrayView_const_iterator(const ArrayView_const_iterator& it)
    : p_(it.p_), arr_(it.arr_), loc_(it.loc_)
{}

template <typename DATA_TYPE, int RANK>
ArrayView_const_iterator<DATA_TYPE,RANK>& ArrayView_const_iterator<DATA_TYPE,RANK>::increment(int d)
{
  ++loc_[d];
  p_ += arr_->stride(d);
  if( loc_[d] == arr_->shape(d) )
  {
    p_ -= arr_->stride(d)*arr_->shape(d);
    loc_[d]=0;
    if(d==0)
    {
      p_ = 0;
      return *this;
    }
    else
    {
      return increment(d-1);
    }
  }
  return *this;
}

//------------------------------------------------------------------------------------------------------

} // namespace atlas

#undef CHECK_RANK
#undef CHECK_BOUNDS
#undef CHECK_BOUNDS_1
#undef CHECK_BOUNDS_2
#undef CHECK_BOUNDS_3
#undef CHECK_BOUNDS_4
#undef CHECK_BOUNDS_5

#endif
