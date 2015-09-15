/*
 * (C) Copyright 1996-2015 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef atlas_ArrayView_h
#error This file can only be included from atlas_ArrayView_h
#endif

#ifndef atlas_ArrayView_impl_h
#define atlas_ArrayView_impl_h

#ifndef ATLAS_ARRAYVIEW_BOUNDS_CHECKING
#include "atlas/atlas_config.h"
#endif

#ifdef ATLAS_ARRAYVIEW_BOUNDS_CHECKING
#include <sstream>
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

//------------------------------------------------------------------------------------------------------

namespace atlas {

//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------

template <typename DATA_TYPE>
ArrayView<DATA_TYPE,0>::ArrayView( const DATA_TYPE* data,
                                   const size_t rank,
                                   const ArrayShape::value_type shape[],
                                   const ArrayStrides::value_type strides[])
  : data_( const_cast<DATA_TYPE*>(data) ), rank_(rank)
{
  strides_.assign(strides,strides+rank_);
  shape_.assign(shape,shape+rank_);
  size_ = 1;
  for( size_t j=0; j<rank; ++j )
    size_ *= shape_[j];
}

template <typename DATA_TYPE>
inline typename ArrayView<DATA_TYPE,0>::iterator ArrayView<DATA_TYPE,0>::begin() { return iterator(this); }

template <typename DATA_TYPE>
inline typename ArrayView<DATA_TYPE,0>::iterator ArrayView<DATA_TYPE,0>::end()   { return iterator(); }

template <typename DATA_TYPE>
inline typename ArrayView<DATA_TYPE,0>::const_iterator ArrayView<DATA_TYPE,0>::cbegin() const { return const_iterator(this); }

template <typename DATA_TYPE>
inline typename ArrayView<DATA_TYPE,0>::const_iterator ArrayView<DATA_TYPE,0>::cend() const { return const_iterator(); }

template <typename DATA_TYPE>
inline typename ArrayView<DATA_TYPE,0>::const_iterator ArrayView<DATA_TYPE,0>::begin() const  { return const_iterator(this); }

template <typename DATA_TYPE>
inline typename ArrayView<DATA_TYPE,0>::const_iterator ArrayView<DATA_TYPE,0>::end()   const  { return const_iterator(); }

template <typename DATA_TYPE>
inline const DATA_TYPE& ArrayView<DATA_TYPE,0>::operator()(size_t i) const { CHECK_RANK(1); CHECK_BOUNDS_1(i); return *(data_+strides_[0]*i); }

template <typename DATA_TYPE>
inline DATA_TYPE&       ArrayView<DATA_TYPE,0>::operator()(size_t i)       { CHECK_RANK(1); CHECK_BOUNDS_1(i); return *(data_+strides_[0]*i); }

template <typename DATA_TYPE>
inline const DATA_TYPE& ArrayView<DATA_TYPE,0>::operator()(size_t i, size_t j) const  { CHECK_RANK(2); CHECK_BOUNDS_2(i,j); return *(data_+strides_[0]*i+j*strides_[1]); }

template <typename DATA_TYPE>
inline DATA_TYPE&       ArrayView<DATA_TYPE,0>::operator()(size_t i, size_t j)        { CHECK_RANK(2); CHECK_BOUNDS_2(i,j); return *(data_+strides_[0]*i+j*strides_[1]); }

template <typename DATA_TYPE>
inline const DATA_TYPE& ArrayView<DATA_TYPE,0>::operator()(size_t i, size_t j, size_t k) const { CHECK_RANK(3); CHECK_BOUNDS_3(i,j,k); return *(data_+strides_[0]*i+j*strides_[1]+k*strides_[2]); }

template <typename DATA_TYPE>
inline DATA_TYPE&       ArrayView<DATA_TYPE,0>::operator()(size_t i, size_t j, size_t k)       { CHECK_RANK(3); CHECK_BOUNDS_3(i,j,k); return *(data_+strides_[0]*i+j*strides_[1]+k*strides_[2]); }

template <typename DATA_TYPE>
inline const DATA_TYPE& ArrayView<DATA_TYPE,0>::operator()(size_t i, size_t j, size_t k, size_t l) const { CHECK_RANK(4); CHECK_BOUNDS_4(i,j,k,l); return *(data_+strides_[0]*i+j*strides_[1]+k*strides_[2]+l*strides_[3]); }

template <typename DATA_TYPE>
inline DATA_TYPE&       ArrayView<DATA_TYPE,0>::operator()(size_t i, size_t j, size_t k, size_t l)       { CHECK_RANK(4); CHECK_BOUNDS_4(i,j,k,l); return *(data_+strides_[0]*i+j*strides_[1]+k*strides_[2]+l*strides_[3]); }

template <typename DATA_TYPE>
inline const DATA_TYPE& ArrayView<DATA_TYPE,0>::operator()(size_t i, size_t j, size_t k, size_t l, size_t m) const { CHECK_RANK(5); CHECK_BOUNDS_5(i,j,k,l,m); return *(data_+strides_[0]*i+j*strides_[1]+k*strides_[2]+l*strides_[3]+m*strides_[4]); }

template <typename DATA_TYPE>
inline DATA_TYPE&       ArrayView<DATA_TYPE,0>::operator()(size_t i, size_t j, size_t k, size_t l, size_t m)       { CHECK_RANK(5); CHECK_BOUNDS_5(i,j,k,l,m); return *(data_+strides_[0]*i+j*strides_[1]+k*strides_[2]+l*strides_[3]+m*strides_[4]); }

template <typename DATA_TYPE>
inline DATA_TYPE& ArrayView<DATA_TYPE,0>::operator()(const ArrayIdx& idx)
{
  CHECK_RANK(idx.size()); CHECK_BOUNDS(idx);
  size_t p=0;
  for( size_t d=0; d<rank(); ++d )
    p += idx[d]*strides_[d];
  return *(data_+p);
}

template <typename DATA_TYPE>
inline const DATA_TYPE& ArrayView<DATA_TYPE,0>::operator()(const ArrayIdx& idx) const
{
  CHECK_RANK(idx.size()); CHECK_BOUNDS(idx);
  size_t p=0;
  for( size_t d=0; d<rank(); ++d )
    p += idx[d]*strides_[d];
  return *(data_+p);
}

template <typename DATA_TYPE>
inline DATA_TYPE*       ArrayView<DATA_TYPE,0>::data()        { return data_; }

template <typename DATA_TYPE>
inline const DATA_TYPE* ArrayView<DATA_TYPE,0>::data() const  { return data_; }

template <typename DATA_TYPE>
inline const ArrayStrides::value_type* ArrayView<DATA_TYPE,0>::strides() const  { return strides_.data(); }

template <typename DATA_TYPE>
inline const ArrayShape::value_type* ArrayView<DATA_TYPE,0>::shape() const  { return shape_.data(); }

template <typename DATA_TYPE>
inline ArrayStrides::value_type ArrayView<DATA_TYPE,0>::stride(size_t i) const { return strides_[i]; }

template <typename DATA_TYPE>
inline ArrayShape::value_type ArrayView<DATA_TYPE,0>::shape(size_t i) const { return shape_[i]; }

template <typename DATA_TYPE>
inline size_t ArrayView<DATA_TYPE,0>::rank() const       { return rank_; }

template <typename DATA_TYPE>
inline size_t ArrayView<DATA_TYPE,0>::size() const       { return size_; }

template <typename DATA_TYPE>
inline void ArrayView<DATA_TYPE,0>::operator=(const DATA_TYPE& scalar) { for(size_t n=0; n<size_; ++n) *(data_+n)=scalar; }


//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------

template <typename DATA_TYPE>
ArrayView<DATA_TYPE,1>::ArrayView( DATA_TYPE* data, const size_t size )
  : data_( const_cast<DATA_TYPE*>(data) )
{
  shape_[0]=size; strides_[0]=1;
}

template <typename DATA_TYPE>
ArrayView<DATA_TYPE,1>::ArrayView( DATA_TYPE* data, const ArrayShape::value_type shape[1], const ArrayStrides::value_type strides[1] )
  : data_( const_cast<DATA_TYPE*>(data) )
{
  strides_[0]=strides[0];       shape_[0]=shape[0];
}

template <typename DATA_TYPE>
ArrayView<DATA_TYPE,1>::ArrayView( const DATA_TYPE* data, const ArrayShape::value_type shape[1] )
  : data_( const_cast<DATA_TYPE*>(data) )
{
  shape_[0]=shape[0]; strides_[0]=1;
}

template <typename DATA_TYPE>
ArrayView<DATA_TYPE,1>::ArrayView( const DATA_TYPE* data, const ArrayShape& shape )
  : data_( const_cast<DATA_TYPE*>(data) )
{
  shape_[0]=shape[0]; strides_[0]=1;
}

template <typename DATA_TYPE>
inline typename ArrayView<DATA_TYPE,1>::iterator ArrayView<DATA_TYPE,1>::begin() { return iterator(this); }

template <typename DATA_TYPE>
inline typename ArrayView<DATA_TYPE,1>::iterator ArrayView<DATA_TYPE,1>::end()   { return iterator(); }

template <typename DATA_TYPE>
inline typename ArrayView<DATA_TYPE,1>::const_iterator ArrayView<DATA_TYPE,1>::cbegin() const { return const_iterator(this); }

template <typename DATA_TYPE>
inline typename ArrayView<DATA_TYPE,1>::const_iterator ArrayView<DATA_TYPE,1>::cend() const  { return const_iterator(); }

template <typename DATA_TYPE>
inline typename ArrayView<DATA_TYPE,1>::const_iterator ArrayView<DATA_TYPE,1>::begin() const  { return const_iterator(this); }

template <typename DATA_TYPE>
inline typename ArrayView<DATA_TYPE,1>::const_iterator ArrayView<DATA_TYPE,1>::end()   const  { return const_iterator(); }

template <typename DATA_TYPE>
inline const DATA_TYPE& ArrayView<DATA_TYPE,1>::operator()(size_t i) const { CHECK_BOUNDS_1(i); return *(data_+strides_[0]*i); }

template <typename DATA_TYPE>
inline DATA_TYPE& ArrayView<DATA_TYPE,1>::operator()(size_t i)       { CHECK_BOUNDS_1(i); return *(data_+strides_[0]*i); }

template <typename DATA_TYPE>
inline const DATA_TYPE& ArrayView<DATA_TYPE,1>::operator()(const ArrayIdx& idx) const
{
  CHECK_RANK(idx.size()); CHECK_BOUNDS(idx);
  size_t p=0;
  for( size_t d=0; d<rank(); ++d )
    p += idx[d]*strides_[d];
  return *(data_+p);
}

template <typename DATA_TYPE>
inline DATA_TYPE& ArrayView<DATA_TYPE,1>::operator()(const ArrayIdx& idx)
{
  CHECK_RANK(idx.size()); CHECK_BOUNDS(idx);
  size_t p=0;
  for( size_t d=0; d<rank(); ++d )
    p += idx[d]*strides_[d];
  return *(data_+p);
}

template <typename DATA_TYPE>
inline const DATA_TYPE& ArrayView<DATA_TYPE,1>::operator[](size_t i) const { CHECK_BOUNDS_1(i); return *(data_+strides_[0]*i); }

template <typename DATA_TYPE>
inline DATA_TYPE& ArrayView<DATA_TYPE,1>::operator[](size_t i)       { CHECK_BOUNDS_1(i); return *(data_+strides_[0]*i); }

template <typename DATA_TYPE>
inline DATA_TYPE* ArrayView<DATA_TYPE,1>::data()        { return data_; }

template <typename DATA_TYPE>
inline const DATA_TYPE* ArrayView<DATA_TYPE,1>::data() const  { return data_; }

template <typename DATA_TYPE>
inline const ArrayStrides::value_type* ArrayView<DATA_TYPE,1>::strides() const   { return strides_; }

template <typename DATA_TYPE>
inline const ArrayShape::value_type* ArrayView<DATA_TYPE,1>::shape() const   { return shape_; }

template <typename DATA_TYPE>
inline ArrayShape::value_type ArrayView<DATA_TYPE,1>::shape(const size_t i) const { return shape_[0]; }

template <typename DATA_TYPE>
inline ArrayStrides::value_type ArrayView<DATA_TYPE,1>::stride(size_t i) const { return strides_[0]; }

template <typename DATA_TYPE>
inline size_t ArrayView<DATA_TYPE,1>::rank() const { return 1; }

template <typename DATA_TYPE>
inline size_t ArrayView<DATA_TYPE,1>::size() const { return shape_[0]; }

template <typename DATA_TYPE>
inline void ArrayView<DATA_TYPE,1>::operator=(const DATA_TYPE& scalar) { for(size_t n=0; n<size(); ++n) *(data_+n)=scalar; }

//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------

template <typename DATA_TYPE>
ArrayView<DATA_TYPE,2>::ArrayView( const DATA_TYPE* data, const ArrayShape::value_type shape[2], const ArrayStrides::value_type strides[2] )
  : data_( const_cast<DATA_TYPE*>(data) )
{
  strides_[0]=strides[0];            shape_[0]=shape[0];
  strides_[1]=strides[1];            shape_[1]=shape[1];
  size_ = shape_[0]*shape_[1];
}

template <typename DATA_TYPE>
ArrayView<DATA_TYPE,2>::ArrayView( const DATA_TYPE* data, const ArrayShape::value_type shape[2] )
  : data_( const_cast<DATA_TYPE*>(data) )
{
  shape_[0]=shape[0]; strides_[0]=shape[1];
  shape_[1]=shape[1]; strides_[1]=1;
  size_ = shape_[0]*shape_[1];
}

template <typename DATA_TYPE>
ArrayView<DATA_TYPE,2>::ArrayView( const DATA_TYPE* data, const std::vector<size_t>& shape )
  : data_( const_cast<DATA_TYPE*>(data) )
{
  shape_[0]=shape[0]; strides_[0]=shape[1];
  shape_[1]=shape[1]; strides_[1]=1;
  size_ = shape_[0]*shape_[1];
}

template <typename DATA_TYPE>
inline typename ArrayView<DATA_TYPE,2>::iterator ArrayView<DATA_TYPE,2>::begin() { return iterator(this); }

template <typename DATA_TYPE>
inline typename ArrayView<DATA_TYPE,2>::iterator ArrayView<DATA_TYPE,2>::end()   { return iterator(); }

template <typename DATA_TYPE>
inline typename ArrayView<DATA_TYPE,2>::const_iterator ArrayView<DATA_TYPE,2>::cbegin() const { return const_iterator(this); }

template <typename DATA_TYPE>
inline typename ArrayView<DATA_TYPE,2>::const_iterator ArrayView<DATA_TYPE,2>::cend() const   { return const_iterator(); }

template <typename DATA_TYPE>
inline typename ArrayView<DATA_TYPE,2>::const_iterator ArrayView<DATA_TYPE,2>::begin() const  { return const_iterator(this); }

template <typename DATA_TYPE>
inline typename ArrayView<DATA_TYPE,2>::const_iterator ArrayView<DATA_TYPE,2>::end()   const  { return const_iterator(); }

template <typename DATA_TYPE>
inline const DATA_TYPE& ArrayView<DATA_TYPE,2>::operator()(size_t i, size_t j) const  { CHECK_BOUNDS_2(i,j); return *(data_+strides_[0]*i+j*strides_[1]); }

template <typename DATA_TYPE>
inline DATA_TYPE&       ArrayView<DATA_TYPE,2>::operator()(size_t i, size_t j)        { CHECK_BOUNDS_2(i,j); return *(data_+strides_[0]*i+j*strides_[1]); }

template <typename DATA_TYPE>
inline const DATA_TYPE& ArrayView<DATA_TYPE,2>::operator()(const ArrayIdx& idx) const
{
  CHECK_RANK(idx.size()); CHECK_BOUNDS(idx);
  size_t p=0;
  for( size_t d=0; d<rank(); ++d )
    p += idx[d]*strides_[d];
  return *(data_+p);
}

template <typename DATA_TYPE>
inline DATA_TYPE& ArrayView<DATA_TYPE,2>::operator()(const ArrayIdx& idx)
{
  CHECK_RANK(idx.size()); CHECK_BOUNDS(idx);
  size_t p=0;
  for( size_t d=0; d<rank(); ++d )
    p += idx[d]*strides_[d];
  return *(data_+p);
}

template <typename DATA_TYPE>
inline const ArrayView<DATA_TYPE,1> ArrayView<DATA_TYPE,2>::operator[](size_t i) const { CHECK_BOUNDS_1(i); return ArrayView<DATA_TYPE,1>( data_+strides_[0]*i, shape_+1, strides_+1 ); }

template <typename DATA_TYPE>
inline ArrayView<DATA_TYPE,1>       ArrayView<DATA_TYPE,2>::operator[](size_t i)       { CHECK_BOUNDS_1(i); return ArrayView<DATA_TYPE,1>( data_+strides_[0]*i, shape_+1, strides_+1 ); }

template <typename DATA_TYPE>
inline DATA_TYPE*       ArrayView<DATA_TYPE,2>::data()            { return data_; }

template <typename DATA_TYPE>
inline const DATA_TYPE* ArrayView<DATA_TYPE,2>::data() const      { return data_; }

template <typename DATA_TYPE>
inline const ArrayStrides::value_type* ArrayView<DATA_TYPE,2>::strides() const   { return strides_; }

template <typename DATA_TYPE>
inline const ArrayShape::value_type* ArrayView<DATA_TYPE,2>::shape() const   { return shape_; }

template <typename DATA_TYPE>
inline ArrayStrides::value_type ArrayView<DATA_TYPE,2>::stride(size_t i) const { return strides_[i]; }

template <typename DATA_TYPE>
inline ArrayShape::value_type ArrayView<DATA_TYPE,2>::shape(size_t i) const {return shape_[i]; }

template <typename DATA_TYPE>
inline size_t ArrayView<DATA_TYPE,2>::rank() const { return 2; }

template <typename DATA_TYPE>
inline size_t ArrayView<DATA_TYPE,2>::size() const { return size_; }

template <typename DATA_TYPE>
inline void ArrayView<DATA_TYPE,2>::operator=(const DATA_TYPE& scalar) { for(size_t n=0; n<size(); ++n) *(data_+n)=scalar; }

//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------

template <typename DATA_TYPE>
ArrayView<DATA_TYPE,3>::ArrayView( const DATA_TYPE* data, const ArrayShape::value_type shape[3], const ArrayStrides::value_type strides [3] ) : data_( const_cast<DATA_TYPE*>(data) )
{
  strides_[0]=strides[0];            shape_[0]=shape[0];
  strides_[1]=strides[1];            shape_[1]=shape[1];
  strides_[2]=strides[2];            shape_[2]=shape[2];
  size_ = shape_[0]*shape_[1]*shape_[2];
}

template <typename DATA_TYPE>
ArrayView<DATA_TYPE,3>::ArrayView( const DATA_TYPE* data, const ArrayShape::value_type shape[3] ) : data_( const_cast<DATA_TYPE*>(data) )
{
  shape_[0]=shape[0]; strides_[0]=shape[2]*shape[1];
  shape_[1]=shape[1]; strides_[1]=shape[2];
  shape_[2]=shape[2]; strides_[2]=1;
  size_ = shape_[0]*shape_[1]*shape_[2];
}

template <typename DATA_TYPE>
ArrayView<DATA_TYPE,3>::ArrayView( const DATA_TYPE* data, const std::vector<size_t>& shape ) : data_( const_cast<DATA_TYPE*>(data) )
{
  shape_[0]=shape[0]; strides_[0]=shape[2]*shape[1];
  shape_[1]=shape[1]; strides_[1]=shape[2];
  shape_[2]=shape[2]; strides_[2]=1;
  size_ = shape_[0]*shape_[1]*shape_[2];
}

template <typename DATA_TYPE>
inline typename ArrayView<DATA_TYPE,3>::iterator ArrayView<DATA_TYPE,3>::begin() { return iterator(this); }

template <typename DATA_TYPE>
inline typename ArrayView<DATA_TYPE,3>::iterator ArrayView<DATA_TYPE,3>::end()   { return iterator(); }

template <typename DATA_TYPE>
inline typename ArrayView<DATA_TYPE,3>::const_iterator ArrayView<DATA_TYPE,3>::cbegin() const { return const_iterator(this); }

template <typename DATA_TYPE>
inline typename ArrayView<DATA_TYPE,3>::const_iterator ArrayView<DATA_TYPE,3>::cend() const   { return const_iterator(); }

template <typename DATA_TYPE>
inline typename ArrayView<DATA_TYPE,3>::const_iterator ArrayView<DATA_TYPE,3>::begin() const  { return const_iterator(this); }

template <typename DATA_TYPE>
inline typename ArrayView<DATA_TYPE,3>::const_iterator ArrayView<DATA_TYPE,3>::end()   const  { return const_iterator(); }

template <typename DATA_TYPE>
inline const DATA_TYPE& ArrayView<DATA_TYPE,3>::operator()(size_t i, size_t j, size_t k) const { CHECK_BOUNDS_3(i,j,k); return *(data_+strides_[0]*i+j*strides_[1]+k*strides_[2]); }

template <typename DATA_TYPE>
inline DATA_TYPE&       ArrayView<DATA_TYPE,3>::operator()(size_t i, size_t j, size_t k)       { CHECK_BOUNDS_3(i,j,k); return *(data_+strides_[0]*i+j*strides_[1]+k*strides_[2]); }

template <typename DATA_TYPE>
inline const DATA_TYPE& ArrayView<DATA_TYPE,3>::operator()(const ArrayIdx& idx) const
{
  CHECK_RANK(idx.size()); CHECK_BOUNDS(idx);
  size_t p=0;
  for( size_t d=0; d<rank(); ++d )
    p += idx[d]*strides_[d];
  return *(data_+p);
}

template <typename DATA_TYPE>
inline DATA_TYPE& ArrayView<DATA_TYPE,3>::operator()(const ArrayIdx& idx)
{
  CHECK_RANK(idx.size()); CHECK_BOUNDS(idx);
  size_t p=0;
  for( size_t d=0; d<rank(); ++d )
    p += idx[d]*strides_[d];
  return *(data_+p);
}

template <typename DATA_TYPE>
inline const ArrayView<DATA_TYPE,2> ArrayView<DATA_TYPE,3>::operator[](size_t i) const { CHECK_BOUNDS_1(i); return ArrayView<DATA_TYPE,2>( data_+strides_[0]*i,shape_+1, strides_+1 ); }

template <typename DATA_TYPE>
inline ArrayView<DATA_TYPE,2>       ArrayView<DATA_TYPE,3>::operator[](size_t i)       { CHECK_BOUNDS_1(i);return ArrayView<DATA_TYPE,2>( data_+strides_[0]*i, shape_+1, strides_+1 ); }

template <typename DATA_TYPE>
inline DATA_TYPE*       ArrayView<DATA_TYPE,3>::data()            { return data_; }

template <typename DATA_TYPE>
inline const DATA_TYPE* ArrayView<DATA_TYPE,3>::data() const      { return data_; }

template <typename DATA_TYPE>
inline const ArrayStrides::value_type* ArrayView<DATA_TYPE,3>::strides() const   { return strides_; }

template <typename DATA_TYPE>
inline const ArrayShape::value_type* ArrayView<DATA_TYPE,3>::shape() const   { return shape_; }

template <typename DATA_TYPE>
inline ArrayStrides::value_type ArrayView<DATA_TYPE,3>::stride(size_t i) const { return strides_[i]; }

template <typename DATA_TYPE>
inline ArrayShape::value_type ArrayView<DATA_TYPE,3>::shape(size_t i) const { return shape_[i]; }

template <typename DATA_TYPE>
inline size_t ArrayView<DATA_TYPE,3>::rank() const { return 3; }

template <typename DATA_TYPE>
inline size_t ArrayView<DATA_TYPE,3>::size() const { return shape_[0]*shape_[1]*shape_[2]; }

template <typename DATA_TYPE>
inline void ArrayView<DATA_TYPE,3>::operator=(const DATA_TYPE& scalar) { for(size_t n=0; n<size(); ++n) *(data_+n)=scalar; }

//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------

template <typename DATA_TYPE>
ArrayView<DATA_TYPE,4>::ArrayView( DATA_TYPE* data, const ArrayShape::value_type shape[4], const ArrayStrides::value_type strides[4] ) : data_( const_cast<DATA_TYPE*>(data) )
{
  strides_[0]=strides[0];            shape_[0]=shape[0];
  strides_[1]=strides[1];            shape_[1]=shape[1];
  strides_[2]=strides[2];            shape_[2]=shape[2];
  strides_[3]=strides[3];            shape_[3]=shape[3];
  size_ = shape_[0]*shape_[1]*shape_[2]*shape_[3];
}

template <typename DATA_TYPE>
ArrayView<DATA_TYPE,4>::ArrayView( const DATA_TYPE* data, const std::vector<size_t>& shape ) : data_( const_cast<DATA_TYPE*>(data) )
{
  shape_[0]=shape[0]; strides_[0]=shape[3]*shape[2]*shape[1];
  shape_[1]=shape[1]; strides_[1]=shape[3]*shape[2];
  shape_[2]=shape[2]; strides_[2]=shape[3];
  shape_[3]=shape[3]; strides_[3]=1;
  size_ = shape_[0]*shape_[1]*shape_[2]*shape_[3];
}

template <typename DATA_TYPE>
ArrayView<DATA_TYPE,4>::ArrayView( const DATA_TYPE* data, const ArrayShape::value_type shape[4] ) : data_( const_cast<DATA_TYPE*>(data) )
{
  shape_[0]=shape[0]; strides_[0]=shape[3]*shape[2]*shape[1];
  shape_[1]=shape[1]; strides_[1]=shape[3]*shape[2];
  shape_[2]=shape[2]; strides_[2]=shape[3];
  shape_[3]=shape[3]; strides_[3]=1;
  size_ = shape_[0]*shape_[1]*shape_[2]*shape_[3];
}

template <typename DATA_TYPE>
inline typename ArrayView<DATA_TYPE,4>::iterator ArrayView<DATA_TYPE,4>::begin() { return iterator(this); }

template <typename DATA_TYPE>
inline typename ArrayView<DATA_TYPE,4>::iterator ArrayView<DATA_TYPE,4>::end()   { return iterator(); }

template <typename DATA_TYPE>
inline typename ArrayView<DATA_TYPE,4>::const_iterator ArrayView<DATA_TYPE,4>::cbegin() const { return const_iterator(this); }

template <typename DATA_TYPE>
inline typename ArrayView<DATA_TYPE,4>::const_iterator ArrayView<DATA_TYPE,4>::cend() const  { return const_iterator(); }

template <typename DATA_TYPE>
inline typename ArrayView<DATA_TYPE,4>::const_iterator ArrayView<DATA_TYPE,4>::begin() const  { return const_iterator(this); }

template <typename DATA_TYPE>
inline typename ArrayView<DATA_TYPE,4>::const_iterator ArrayView<DATA_TYPE,4>::end()   const  { return const_iterator(); }

template <typename DATA_TYPE>
inline const DATA_TYPE& ArrayView<DATA_TYPE,4>::operator()(size_t i, size_t j, size_t k, size_t l) const { CHECK_BOUNDS_4(i,j,k,l); return *(data_+strides_[0]*i+j*strides_[1]+k*strides_[2]+l*strides_[3]); }

template <typename DATA_TYPE>
inline DATA_TYPE&       ArrayView<DATA_TYPE,4>::operator()(size_t i, size_t j, size_t k, size_t l)       { CHECK_BOUNDS_4(i,j,k,l); return *(data_+strides_[0]*i+j*strides_[1]+k*strides_[2]+l*strides_[3]); }

template <typename DATA_TYPE>
inline const DATA_TYPE& ArrayView<DATA_TYPE,4>::operator()(const ArrayIdx& idx) const
{
  CHECK_RANK(idx.size()); CHECK_BOUNDS(idx);
  size_t p=0;
  for( size_t d=0; d<rank(); ++d )
    p += idx[d]*strides_[d];
  return *(data_+p);
}

template <typename DATA_TYPE>
inline DATA_TYPE& ArrayView<DATA_TYPE,4>::operator()(const ArrayIdx& idx)
{
  CHECK_RANK(idx.size()); CHECK_BOUNDS(idx);
  size_t p=0;
  for( size_t d=0; d<rank(); ++d )
    p += idx[d]*strides_[d];
  return *(data_+p);
}

template <typename DATA_TYPE>
inline const ArrayView<DATA_TYPE,3> ArrayView<DATA_TYPE,4>::operator[](size_t i) const { CHECK_BOUNDS_1(i); return ArrayView<DATA_TYPE,3>( data_+strides_[0]*i, shape_+1, strides_+1 ); }

template <typename DATA_TYPE>
inline ArrayView<DATA_TYPE,3>       ArrayView<DATA_TYPE,4>::operator[](size_t i)       { CHECK_BOUNDS_1(i); return ArrayView<DATA_TYPE,3>( data_+strides_[0]*i, shape_+1, strides_+1 ); }

template <typename DATA_TYPE>
inline DATA_TYPE*       ArrayView<DATA_TYPE,4>::data()            { return data_; }

template <typename DATA_TYPE>
inline const DATA_TYPE* ArrayView<DATA_TYPE,4>::data() const      { return data_; }

template <typename DATA_TYPE>
inline const ArrayStrides::value_type* ArrayView<DATA_TYPE,4>::strides() const   { return strides_; }

template <typename DATA_TYPE>
inline const ArrayShape::value_type* ArrayView<DATA_TYPE,4>::shape() const   { return shape_; }

template <typename DATA_TYPE>
inline ArrayStrides::value_type ArrayView<DATA_TYPE,4>::stride(size_t i) const { return strides_[i]; }

template <typename DATA_TYPE>
inline ArrayShape::value_type ArrayView<DATA_TYPE,4>::shape(size_t i) const { return shape_[i]; }

template <typename DATA_TYPE>
inline size_t ArrayView<DATA_TYPE,4>::rank() const { return 4; }

template <typename DATA_TYPE>
inline size_t ArrayView<DATA_TYPE,4>::size() const { return shape_[0]*shape_[1]*shape_[2]*shape_[3]; }

template <typename DATA_TYPE>
inline void ArrayView<DATA_TYPE,4>::operator=(const DATA_TYPE& scalar) { for(size_t n=0; n<size(); ++n) *(data_+n)=scalar; }


//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------
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
//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------
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
