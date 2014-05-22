/*
 * (C) Copyright 1996-2014 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */



/// @file ArrayView.hpp
/// This file contains the ArrayView class, a class that allows to wrap any contiguous raw data into
/// a view which is accessible with multiple indices.
/// All it needs is the strides for each index, and the extent of each index.
/// ATTENTION: The last index is stride 1
///
/// Bounds-checking can be turned ON by defining "ATLAS_ARRAYVIEW_BOUNDS_CHECKING"
/// before including this header.
///  
/// Example:
/// int[] array = { 1, 2, 3, 4, 5, 6, 7, 8, 9};
/// int[2] strides = { 3, 1 };
/// int[2] extents = { 3, 3 };
/// ArrayView<int,2> matrix( array, strides, extents );
/// for( int i=0; i<matrix.extent(0); ++i ) {
///   for( int j=0; j<matrix.extent(1); ++j ) {
///     matrix(i,j) *= 10;
///   }
/// }
///
/// There is also an easier way to wrap Field and Array classes:
/// ArrayView<int,3> fieldview( Field );
/// ArrayView<int,2> arrayview( Array );

#ifndef atlas_ArrayView_hpp
#define atlas_ArrayView_hpp

// #define ATLAS_ARRAYVIEW_BOUNDS_CHECKING

#ifdef ATLAS_ARRAYVIEW_BOUNDS_CHECKING
#include <stdexcept>
#define CHECK_BOUNDS_1(i)\
  if(i>=extents_[0]) {throw std::range_error("index 'i' out of bounds");}
#define CHECK_BOUNDS_2(i,j)\
  if(i>=extents_[0]) {throw std::range_error("index 'i' out of bounds");} \
  if(j>=extents_[1]) {throw std::range_error("index 'j' out of bounds");}
#define CHECK_BOUNDS_3(i,j,k)\
  if(i>=extents_[0]) {throw std::range_error("index 'i' out of bounds");} \
  if(j>=extents_[1]) {throw std::range_error("index 'j' out of bounds");} \
  if(k>=extents_[2]) {throw std::range_error("index 'k' out of bounds");}
#define CHECK_BOUNDS_4(i,j,k,l)\
  if(i>=extents_[0]) {throw std::range_error("index 'i' out of bounds");} \
  if(j>=extents_[1]) {throw std::range_error("index 'j' out of bounds");} \
  if(k>=extents_[2]) {throw std::range_error("index 'k' out of bounds");} \
  if(l>=extents_[3]) {throw std::range_error("index 'l' out of bounds");}  
#else
#define CHECK_BOUNDS_1(i)
#define CHECK_BOUNDS_2(i,j)
#define CHECK_BOUNDS_3(i,j,k)
#define CHECK_BOUNDS_4(i,j,k,l)
#endif


//------------------------------------------------------------------------------------------------------

namespace atlas {

  class Field;
  template< typename DATA_TYPE > class Array;

//------------------------------------------------------------------------------------------------------

template< typename DATA_TYPE, int RANK >
class ArrayView
{};

//------------------------------------------------------------------------------------------------------

template< typename DATA_TYPE >
class ArrayView < DATA_TYPE, 1 >
{
public:
  ArrayView() {}
  ArrayView( DATA_TYPE* data, const int strides[1], const int extents[1] ) : data_( const_cast<DATA_TYPE*>(data) )
  {
    strides_[0]=strides[0];       extents_[0]=extents[0];
  }
  ArrayView( const Array<DATA_TYPE>& array );
  ArrayView( const Field& field );
  
  const DATA_TYPE& operator()(int i) const { CHECK_BOUNDS_1(i); return *(data_+strides_[0]*i); }
  DATA_TYPE&       operator()(int i)       { CHECK_BOUNDS_1(i); return *(data_+strides_[0]*i); }
  const DATA_TYPE& operator[](int i) const { CHECK_BOUNDS_1(i); return *(data_+strides_[0]*i); }
  DATA_TYPE&       operator[](int i)       { CHECK_BOUNDS_1(i); return *(data_+strides_[0]*i); }
  DATA_TYPE* data()            { return data_; }
  const int* strides() const   { return strides_; }
  const int* extents() const   { return extents_; }
  std::size_t size() const { return extents_[0]; }
  std::size_t total_size() const { return extents_[0]; }
  void operator=(const DATA_TYPE& scalar) { for(int n=0; n<total_size(); ++n) *(data_+n)=scalar; }
private:
  DATA_TYPE* data_;
  int strides_[1];
  int extents_[1];
};

//------------------------------------------------------------------------------------------------------

template< typename DATA_TYPE >
class ArrayView < DATA_TYPE, 2 >
{
public:
  
  ArrayView() {}
  ArrayView( const DATA_TYPE* data, const int strides[2], const int extents[2] ) : data_( const_cast<DATA_TYPE*>(data) )
  {
    strides_[0]=strides[0];            extents_[0]=extents[0];
    strides_[1]=strides[1];            extents_[1]=extents[1];
  }
  ArrayView( const Array<DATA_TYPE>& array );
  ArrayView( const Field& field );
  
  const DATA_TYPE& operator()(int i, int j) const  { CHECK_BOUNDS_2(i,j); return *(data_+strides_[0]*i+j*strides_[1]); }
  DATA_TYPE&       operator()(int i, int j)        { CHECK_BOUNDS_2(i,j); return *(data_+strides_[0]*i+j*strides_[1]); }
  const ArrayView<DATA_TYPE,1> operator[](int i) const { CHECK_BOUNDS_1(i);return ArrayView<DATA_TYPE,1>( data_+strides_[0]*i, strides_+1, extents_+1 ); }
  ArrayView<DATA_TYPE,1>       operator[](int i)       { return ArrayView<DATA_TYPE,1>( data_+strides_[0]*i, strides_+1, extents_+1 ); }
  DATA_TYPE* data()            { return data_; }
  const int* strides() const   { return strides_; }
  const int* extents() const   { return extents_; }
  std::size_t size() const { return extents_[0]; }
  std::size_t total_size() const { return extents_[0]*extents_[1]; }
  void operator=(const DATA_TYPE& scalar) { for(int n=0; n<total_size(); ++n) *(data_+n)=scalar; }
private:
  DATA_TYPE* data_;
  int strides_[2];
  int extents_[2];
};

//------------------------------------------------------------------------------------------------------

template< typename DATA_TYPE >
class ArrayView < DATA_TYPE, 3 >
{
public:
  ArrayView() {}
  ArrayView( const DATA_TYPE* data, const int strides[3], const int extents[3] ) : data_( const_cast<DATA_TYPE*>(data) )
  {
    strides_[0]=strides[0];            extents_[0]=extents[0];
    strides_[1]=strides[1];            extents_[1]=extents[1];
    strides_[2]=strides[2];            extents_[2]=extents[2];
  }
  ArrayView( const Array<DATA_TYPE>& array );
  ArrayView( const Field& field );
  
  const DATA_TYPE& operator()(int i, int j, int k) const { CHECK_BOUNDS_3(i,j,k); return *(data_+strides_[0]*i+j*strides_[1]+k*strides_[2]); }
  DATA_TYPE&       operator()(int i, int j, int k)       { CHECK_BOUNDS_3(i,j,k); return *(data_+strides_[0]*i+j*strides_[1]+k*strides_[2]); }
  const ArrayView<DATA_TYPE,2> operator[](int i) const { CHECK_BOUNDS_1(i); return ArrayView<DATA_TYPE,2>( data_+strides_[0]*i, strides_+1, extents_+1 ); }
  ArrayView<DATA_TYPE,2>       operator[](int i)       { CHECK_BOUNDS_1(i);return ArrayView<DATA_TYPE,2>( data_+strides_[0]*i, strides_+1, extents_+1 ); }
  DATA_TYPE* data()            { return data_; }
  const int* strides() const   { return strides_; }
  const int* extents() const   { return extents_; }
  std::size_t size() const { return extents_[0]; }
  std::size_t total_size() const { return extents_[0]*extents_[1]*extents_[2]; }
  void operator=(const DATA_TYPE& scalar) { for(int n=0; n<total_size(); ++n) *(data_+n)=scalar; }
private:
  DATA_TYPE* data_;
  int strides_[3];
  int extents_[3];
};

//------------------------------------------------------------------------------------------------------

template< typename DATA_TYPE >
class ArrayView < DATA_TYPE, 4 >
{
public:
  ArrayView() {}
  ArrayView( DATA_TYPE* data, const int strides[4], const int extents[4] ) : data_( const_cast<DATA_TYPE*>(data) )
  {
    strides_[0]=strides[0];            extents_[0]=extents[0];
    strides_[1]=strides[1];            extents_[1]=extents[1];
    strides_[2]=strides[2];            extents_[2]=extents[2];
    strides_[3]=strides[3];            extents_[3]=extents[3];
  }
  ArrayView( const Array<DATA_TYPE>& array );
  ArrayView( const Field& field );
  
  const DATA_TYPE& operator()(int i, int j, int k, int l) const { CHECK_BOUNDS_4(i,j,k,l); return *(data_+strides_[0]*i+j*strides_[1]+k*strides_[2]+l*strides_[3]); }
  DATA_TYPE&       operator()(int i, int j, int k, int l)       { CHECK_BOUNDS_4(i,j,k,l); return *(data_+strides_[0]*i+j*strides_[1]+k*strides_[2]+l*strides_[3]); }
  const ArrayView<DATA_TYPE,3> operator[](int i) const { CHECK_BOUNDS_1(i); return ArrayView<DATA_TYPE,3>( data_+strides_[0]*i, strides_+1, extents_+1 ); }
  ArrayView<DATA_TYPE,3>       operator[](int i)       { CHECK_BOUNDS_1(i); return ArrayView<DATA_TYPE,3>( data_+strides_[0]*i, strides_+1, extents_+1 ); }
  DATA_TYPE* data()            { return data_; }
  const int* strides() const   { return strides_; }
  const int* extents() const   { return extents_; }
  std::size_t size() const { return extents_[0]; }
  std::size_t total_size() const { return extents_[0]*extents_[1]*extents_[2]*extents_[3]; }
  void operator=(const DATA_TYPE& scalar) { for(int n=0; n<total_size(); ++n) *(data_+n)=scalar; }
  
private:
  DATA_TYPE* data_;
  int strides_[4];
  int extents_[4];
};

//------------------------------------------------------------------------------------------------------

} // namespace atlas

#undef CHECK_BOUNDS_1
#undef CHECK_BOUNDS_2
#undef CHECK_BOUNDS_3
#undef CHECK_BOUNDS_4

#endif
