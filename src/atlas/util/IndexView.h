/*
 * (C) Copyright 1996-2015 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */



/// @file IndexView.h
/// This file contains the IndexView class, a class that allows to wrap any contiguous raw data into
/// a view which is accessible with multiple indices.
/// This view is intended to work with Connectivity Tables storing Fortran Numbering internally
/// All it needs is the strides for each index, and the shape of each index.
/// ATTENTION: The last index is stride 1
///
/// Bounds-checking can be turned ON by defining "ATLAS_INDEXVIEW_BOUNDS_CHECKING"
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

#ifndef atlas_IndexView_h
#define atlas_IndexView_h

#include "atlas/atlas_defines.h"

#define ATLAS_INDEXVIEW_BOUNDS_CHECKING

#ifdef ATLAS_INDEXVIEW_BOUNDS_CHECKING
#include <eckit/exception/Exceptions.h>


#define CHECK_RANK(R)\
  if(rank()!=R) { std::ostringstream msg; msg << "IndexView  rank mismatch: rank()="<<rank()<< " != " << R; throw eckit::OutOfRange(msg.str(),Here()); }
#define CHECK_BOUNDS(idx) {\
  for( size_t d=0; d<rank(); ++d ) { \
    if(idx[d]>=shape_[d]) {std::ostringstream msg; msg << "index " << d << " out of bounds: " << idx[d] << " >= " << shape_[d]; throw eckit::OutOfRange(msg.str(),Here()); } } }
#define CHECK_BOUNDS_1(i)\
	if(i>=shape_[0]) {std::ostringstream msg; msg << "IndexView(i) index out of bounds: i=" << i << " >= " << shape_[0]; throw eckit::OutOfRange(msg.str(),Here()); }
#define CHECK_BOUNDS_2(i,j)\
	if(i>=shape_[0]) {std::ostringstream msg; msg << "IndexView(i,j) index out of bounds: i=" << i << " >= " << shape_[0]; throw eckit::OutOfRange(msg.str(),Here()); }\
	if(j>=shape_[1]) {std::ostringstream msg; msg << "IndexView(i,j) index out of bounds: j=" << j << " >= " << shape_[1]; throw eckit::OutOfRange(msg.str(),Here()); }
#define CHECK_BOUNDS_3(i,j,k)\
	if(i>=shape_[0]) {std::ostringstream msg; msg << "IndexView(i,j,k) index out of bounds: i=" << i << " >= " << shape_[0]; throw eckit::OutOfRange(msg.str(),Here()); }\
	if(j>=shape_[1]) {std::ostringstream msg; msg << "IndexView(i,j,k) index out of bounds: j=" << j << " >= " << shape_[1]; throw eckit::OutOfRange(msg.str(),Here()); }\
	if(k>=shape_[2]) {std::ostringstream msg; msg << "IndexView(i,j,k) index out of bounds: k=" << k << " >= " << shape_[2]; throw eckit::OutOfRange(msg.str(),Here()); }
#define CHECK_BOUNDS_4(i,j,k,l)\
	if(i>=shape_[0]) {std::ostringstream msg; msg << "IndexView(i,j,k,l) index out of bounds: i=" << i << " >= " << shape_[0]; throw eckit::OutOfRange(msg.str(),Here()); }\
	if(j>=shape_[1]) {std::ostringstream msg; msg << "IndexView(i,j,k,l) index out of bounds: j=" << j << " >= " << shape_[1]; throw eckit::OutOfRange(msg.str(),Here()); }\
	if(k>=shape_[2]) {std::ostringstream msg; msg << "IndexView(i,j,k,l) index out of bounds: k=" << k << " >= " << shape_[2]; throw eckit::OutOfRange(msg.str(),Here()); }\
	if(l>=shape_[3]) {std::ostringstream msg; msg << "IndexView(i,j,k,l) index out of bounds: l=" << l << " >= " << shape_[3]; throw eckit::OutOfRange(msg.str(),Here()); }
#define CHECK_BOUNDS_5(i,j,k,l,m)\
	if(i>=shape_[0]) {std::ostringstream msg; msg << "IndexView(i,j,k,l,m) index out of bounds: i=" << i << " >= " << shape_[0]; throw eckit::OutOfRange(msg.str(),Here()); }\
	if(j>=shape_[1]) {std::ostringstream msg; msg << "IndexView(i,j,k,l,m) index out of bounds: j=" << j << " >= " << shape_[1]; throw eckit::OutOfRange(msg.str(),Here()); }\
	if(k>=shape_[2]) {std::ostringstream msg; msg << "IndexView(i,j,k,l,m) index out of bounds: k=" << k << " >= " << shape_[2]; throw eckit::OutOfRange(msg.str(),Here()); }\
	if(l>=shape_[3]) {std::ostringstream msg; msg << "IndexView(i,j,k,l,m) index out of bounds: l=" << l << " >= " << shape_[3]; throw eckit::OutOfRange(msg.str(),Here()); }\
	if(m>=shape_[4]) {std::ostringstream msg; msg << "IndexView(i,j,k,l,m) index out of bounds: m=" << m << " >= " << shape_[4]; throw eckit::OutOfRange(msg.str(),Here()); }
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

  class Field;
  template< typename DATA_TYPE > class ArrayT;

//------------------------------------------------------------------------------------------------------

template< typename DATA_TYPE, int RANK >
class IndexView
{};

//------------------------------------------------------------------------------------------------------

namespace detail {
// FortranIndex:
// Helper class that does +1 and -1 operations on stored values
template< typename DATA_TYPE >
class FortranIndex
{
public:
  enum { BASE = 1 };
public:
  FortranIndex(DATA_TYPE* idx): idx_(idx) {}
  void set(const DATA_TYPE& value) { *(idx_) = value+BASE; }
  DATA_TYPE get() const { return *(idx_)-BASE; }
  void operator=(const DATA_TYPE& value) { set(value); }
  FortranIndex<DATA_TYPE>& operator=(const FortranIndex<DATA_TYPE>& other) { set(other.get()); return *this; }
  FortranIndex<DATA_TYPE>& operator+(const DATA_TYPE& value) { *(idx_)+=value; return *this; }
  FortranIndex<DATA_TYPE>& operator-(const DATA_TYPE& value) { *(idx_)-=value; return *this; }
  FortranIndex<DATA_TYPE>& operator--() { --(*(idx_)); return *this; }
  FortranIndex<DATA_TYPE>& operator++() { ++(*(idx_)); return *this; }

  //implicit conversion
  operator DATA_TYPE() const { return get(); }

private:
  DATA_TYPE* idx_;
};
}


#ifdef ATLAS_HAVE_FORTRAN
#define INDEX_REF Index
#define FROM_FORTRAN -1
#define TO_FORTRAN   +1
#else
#define INDEX_REF *
#define FROM_FORTRAN
#define TO_FORTRAN
#endif

//------------------------------------------------------------------------------------------------------

template< typename DATA_TYPE >
class IndexView < DATA_TYPE, 1 >
{
public:
#ifdef ATLAS_HAVE_FORTRAN
  typedef detail::FortranIndex<DATA_TYPE> Index;
#else
    typedef DATA_TYPE& Index;
#endif

public:
  IndexView() {}
  IndexView( DATA_TYPE* data, const size_t strides[1], const size_t shape[1] ) : data_( const_cast<DATA_TYPE*>(data) )
  {
    strides_[0]=strides[0];       shape_[0]=shape[0];
  }
  IndexView( const ArrayT<DATA_TYPE>& array );
  IndexView( const Field& field );

  DATA_TYPE operator()(size_t i) const { CHECK_BOUNDS_1(i); return *(data_+strides_[0]*i) FROM_FORTRAN; }
  Index     operator()(size_t i)       { CHECK_BOUNDS_1(i); return INDEX_REF(data_+strides_[0]*i); }

  DATA_TYPE operator[](size_t i) const { CHECK_BOUNDS_1(i); return *(data_+strides_[0]*i) FROM_FORTRAN; }
  Index     operator[](size_t i)       { CHECK_BOUNDS_1(i); return INDEX_REF(data_+strides_[0]*i); }

  const size_t* strides() const   { return strides_; }
  const size_t* shape() const   { return shape_; }
  size_t shape(const size_t i) const { return shape_[0]; }

  size_t size() const { return shape_[0]; }
  size_t total_size() const { return shape_[0]; }
  void operator=(const DATA_TYPE& scalar) { for(size_t n=0; n<total_size(); ++n) *(data_+n)=scalar TO_FORTRAN; }
private:
  DATA_TYPE* data_;
  size_t strides_[1];
  size_t shape_[1];
};

//------------------------------------------------------------------------------------------------------

template< typename DATA_TYPE >
class IndexView < DATA_TYPE, 2 >
{
public:
#ifdef ATLAS_HAVE_FORTRAN
    typedef detail::FortranIndex<DATA_TYPE> Index;
#else
    typedef DATA_TYPE& Index;
#endif

public:

  IndexView() {}
  IndexView( const DATA_TYPE* data, const size_t strides[2], const size_t shape[2] ) : data_( const_cast<DATA_TYPE*>(data) )
  {
    strides_[0]=strides[0];            shape_[0]=shape[0];
    strides_[1]=strides[1];            shape_[1]=shape[1];
  }
  IndexView( const DATA_TYPE* data, const size_t shape[2] ) : data_( const_cast<DATA_TYPE*>(data) )
  {
    shape_[0]=shape[0]; strides_[0]=shape[1];
    shape_[1]=shape[1]; strides_[1]=1;
  }
  IndexView( const ArrayT<DATA_TYPE>& array );
  IndexView( const Field& field );

  DATA_TYPE operator()(size_t i, size_t j) const  { CHECK_BOUNDS_2(i,j); return *(data_+strides_[0]*i+j*strides_[1]) FROM_FORTRAN; }
  Index     operator()(size_t i, size_t j)        { CHECK_BOUNDS_2(i,j); return INDEX_REF(data_+strides_[0]*i+j*strides_[1]); }


  const IndexView<DATA_TYPE,1> operator[](size_t i) const { CHECK_BOUNDS_1(i); return IndexView<DATA_TYPE,1>( data_+strides_[0]*i, strides_+1, shape_+1 ); }
  IndexView<DATA_TYPE,1>       operator[](size_t i)       { CHECK_BOUNDS_1(i); return IndexView<DATA_TYPE,1>( data_+strides_[0]*i, strides_+1, shape_+1 ); }

  const size_t* strides() const   { return strides_; }
  const size_t* shape() const   { return shape_; }
  size_t shape(const size_t i) const { return shape_[i]; }

  size_t size() const { return shape_[0]; }
  size_t total_size() const { return shape_[0]*shape_[1]; }
  void operator=(const DATA_TYPE& scalar) { for(size_t n=0; n<total_size(); ++n) *(data_+n)=scalar TO_FORTRAN; }
private:
  DATA_TYPE* data_;
  size_t strides_[2];
  size_t shape_[2];
};

//------------------------------------------------------------------------------------------------------

template< typename DATA_TYPE >
class IndexView < DATA_TYPE, 3 >
{
public:
#ifdef ATLAS_HAVE_FORTRAN
    typedef detail::FortranIndex<DATA_TYPE> Index;
#else
    typedef DATA_TYPE& Index;
#endif

public:
  IndexView() {}
  IndexView( const DATA_TYPE* data, const size_t strides[3], const size_t shape[3] ) : data_( const_cast<DATA_TYPE*>(data) )
  {
    strides_[0]=strides[0];            shape_[0]=shape[0];
    strides_[1]=strides[1];            shape_[1]=shape[1];
    strides_[2]=strides[2];            shape_[2]=shape[2];
  }
  IndexView( const DATA_TYPE* data, const size_t shape[3] ) : data_( const_cast<DATA_TYPE*>(data) )
  {
    shape_[0]=shape[0]; strides_[0]=shape[2]*shape[1];
    shape_[1]=shape[1]; strides_[1]=shape[2];
    shape_[2]=shape[2]; strides_[2]=1;
  }
  IndexView( const DATA_TYPE* data, const std::vector<size_t>& shape ) : data_( const_cast<DATA_TYPE*>(data) )
  {
    shape_[0]=shape[0]; strides_[0]=shape[2]*shape[1];
    shape_[1]=shape[1]; strides_[1]=shape[2];
    shape_[2]=shape[2]; strides_[2]=1;
  }
  IndexView( const ArrayT<DATA_TYPE>& array );
  IndexView( const Field& field );

  DATA_TYPE operator()(size_t i, size_t j, size_t k) const { CHECK_BOUNDS_3(i,j,k); return *(data_+strides_[0]*i+j*strides_[1]+k*strides_[2]) FROM_FORTRAN; }
  Index     operator()(size_t i, size_t j, size_t k)       { CHECK_BOUNDS_3(i,j,k); return INDEX_REF(data_+strides_[0]*i+j*strides_[1]+k*strides_[2]); }

  const IndexView<DATA_TYPE,2> operator[](size_t i) const { CHECK_BOUNDS_1(i); return IndexView<DATA_TYPE,2>( data_+strides_[0]*i, strides_+1, shape_+1 ); }
  IndexView<DATA_TYPE,2>       operator[](size_t i)       { CHECK_BOUNDS_1(i); return IndexView<DATA_TYPE,2>( data_+strides_[0]*i, strides_+1, shape_+1 ); }

  const size_t* strides() const   { return strides_; }
  const size_t* shape() const   { return shape_; }
  size_t shape(const size_t i) const { return shape_[i]; }

  size_t size() const { return shape_[0]; }
  size_t total_size() const { return shape_[0]*shape_[1]*shape_[2]; }
  void operator=(const DATA_TYPE& scalar) { for(size_t n=0; n<total_size(); ++n) *(data_+n)=scalar TO_FORTRAN; }
private:
  DATA_TYPE* data_;
  size_t strides_[3];
  size_t shape_[3];
};

//------------------------------------------------------------------------------------------------------

template< typename DATA_TYPE >
class IndexView < DATA_TYPE, 4 >
{
public:
#ifdef ATLAS_HAVE_FORTRAN
    typedef detail::FortranIndex<DATA_TYPE> Index;
#else
    typedef DATA_TYPE& Index;
#endif

public:
  IndexView() {}
  IndexView( DATA_TYPE* data, const size_t strides[4], const size_t shape[4] ) : data_( const_cast<DATA_TYPE*>(data) )
  {
    strides_[0]=strides[0];            shape_[0]=shape[0];
    strides_[1]=strides[1];            shape_[1]=shape[1];
    strides_[2]=strides[2];            shape_[2]=shape[2];
    strides_[3]=strides[3];            shape_[3]=shape[3];
  }
  IndexView( const ArrayT<DATA_TYPE>& array );
  IndexView( const Field& field );

  DATA_TYPE operator()(size_t i, size_t j, size_t k, size_t l) const { CHECK_BOUNDS_4(i,j,k,l); return *(data_+strides_[0]*i+j*strides_[1]+k*strides_[2]+l*strides_[3]) FROM_FORTRAN; }
  Index     operator()(size_t i, size_t j, size_t k, size_t l)       { CHECK_BOUNDS_4(i,j,k,l); return INDEX_REF(data_+strides_[0]*i+j*strides_[1]+k*strides_[2]+l*strides_[3]); }

  const IndexView<DATA_TYPE,3> operator[](size_t i) const { CHECK_BOUNDS_1(i); return IndexView<DATA_TYPE,3>( data_+strides_[0]*i, strides_+1, shape_+1 ); }
  IndexView<DATA_TYPE,3>       operator[](size_t i)       { CHECK_BOUNDS_1(i); return IndexView<DATA_TYPE,3>( data_+strides_[0]*i, strides_+1, shape_+1 ); }

  const size_t* strides() const   { return strides_; }
  const size_t* shape() const   { return shape_; }
  size_t shape(const size_t i) const { return shape_[i]; }

  size_t size() const { return shape_[0]; }
  size_t total_size() const { return shape_[0]*shape_[1]*shape_[2]*shape_[3]; }
  void operator=(const DATA_TYPE& scalar) { for(size_t n=0; n<total_size(); ++n) *(data_+n)=scalar TO_FORTRAN; }

private:
  DATA_TYPE* data_;
  size_t strides_[4];
  size_t shape_[4];
};

//------------------------------------------------------------------------------------------------------

} // namespace atlas

#undef CHECK_RANK
#undef CHECK_BOUNDS
#undef CHECK_BOUNDS_1
#undef CHECK_BOUNDS_2
#undef CHECK_BOUNDS_3
#undef CHECK_BOUNDS_4
#undef CHECK_BOUNDS_5
#undef FROM_FORTRAN
#undef TO_FORTRAN
#undef INDEX_REF
#endif
