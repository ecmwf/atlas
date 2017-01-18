/*
 * (C) Copyright 1996-2016 ECMWF.
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
/// IndexView<int,3> fieldview( field::Field );
/// IndexView<int,2> INDEXVIEW( Array );

#ifndef atlas_IndexView_h
#define atlas_IndexView_h

#include "atlas/internals/atlas_config.h"
#include "atlas/internals/atlas_defines.h"
#ifdef ATLAS_HAVE_GRIDTOOLS_STORAGE
#include "atlas/array/gridtools/GridToolsIndexView.h"
#else
#include <iosfwd>

#ifdef ATLAS_INDEXVIEW_BOUNDS_CHECKING
#include <eckit/exception/Exceptions.h>


#define CHECK_RANK(R)\
  if(Rank!=R) { std::ostringstream msg; msg << "IndexView  rank mismatch: Rank="<<Rank<< " != " << R; throw eckit::OutOfRange(msg.str(),Here()); }
#define CHECK_BOUNDS_1(i)\
  if(i>=shape_[0]) {std::ostringstream msg; msg << "IndexView(i) index out of bounds: i=" << i << " >= " << shape_[0]; throw eckit::OutOfRange(msg.str(),Here()); }
#define CHECK_BOUNDS_2(i,j)\
	if(i>=shape_[0]) {std::ostringstream msg; msg << "IndexView(i,j) index out of bounds: i=" << i << " >= " << shape_[0]; throw eckit::OutOfRange(msg.str(),Here()); }\
	if(j>=shape_[1]) {std::ostringstream msg; msg << "IndexView(i,j) index out of bounds: j=" << j << " >= " << shape_[1]; throw eckit::OutOfRange(msg.str(),Here()); }
#else
#define CHECK_RANK(R)
#define CHECK_BOUNDS_1(i)
#define CHECK_BOUNDS_2(i,j)
#endif

//------------------------------------------------------------------------------------------------------

namespace atlas {
namespace array {

//------------------------------------------------------------------------------------------------------

namespace detail {
// FortranIndex:
// Helper class that does +1 and -1 operations on stored values
template< typename Value >
class FortranIndex
{
public:
  enum { BASE = 1 };
public:
  FortranIndex(Value* idx): idx_(idx) {}
  void set(const Value& value) { *(idx_) = value+BASE; }
  Value get() const { return *(idx_)-BASE; }
  void operator=(const Value& value) { set(value); }
  FortranIndex<Value>& operator=(const FortranIndex<Value>& other) { set(other.get()); return *this; }
  FortranIndex<Value>& operator+(const Value& value) { *(idx_)+=value; return *this; }
  FortranIndex<Value>& operator-(const Value& value) { *(idx_)-=value; return *this; }
  FortranIndex<Value>& operator--() { --(*(idx_)); return *this; }
  FortranIndex<Value>& operator++() { ++(*(idx_)); return *this; }

  //implicit conversion
  operator Value() const { return get(); }

private:
  Value* idx_;
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

template< typename Value, int Rank >
class IndexView {
public:
#ifdef ATLAS_HAVE_FORTRAN
  typedef detail::FortranIndex<Value> Index;
#else
  typedef Value& Index;
#endif

public:

  IndexView( Value* data, const size_t shape[Rank] );

  Value operator()(size_t i) const { CHECK_RANK(1); CHECK_BOUNDS_1(i); return *(data_+strides_[0]*i) FROM_FORTRAN; }
  Index operator()(size_t i)       { CHECK_RANK(1); CHECK_BOUNDS_1(i); return INDEX_REF(data_+strides_[0]*i); }
  Value operator()(size_t i, size_t j) const { CHECK_RANK(2); CHECK_BOUNDS_2(i,j); return *(data_+strides_[0]*i+strides_[1]*j) FROM_FORTRAN; }
  Index operator()(size_t i, size_t j)       { CHECK_RANK(2); CHECK_BOUNDS_2(i,j); return INDEX_REF(data_+strides_[0]*i+strides_[1]*j); }

  size_t size() const { return shape_[0]; }

  void dump(std::ostream& os) const;

private:
  Value* data_;
  size_t strides_[Rank];
  size_t shape_[Rank];
};

//------------------------------------------------------------------------------------------------------

} // namespace array
} // namespace atlas

#endif

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
