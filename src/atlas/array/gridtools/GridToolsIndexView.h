/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#pragma once

#include "atlas/internals/atlas_defines.h"
#include "atlas/array/gridtools/GridToolsTraits.h"

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

//------------------------------------------------------------------------------------------------------

template< typename Value, int Rank >
class IndexView
{
public:
// -- Type definitions
#ifdef ATLAS_HAVE_FORTRAN
    typedef detail::FortranIndex<Value> Index;
    #define INDEX_REF Index
    #define FROM_FORTRAN -1
    #define TO_FORTRAN   +1
#else
    typedef Value& Index;
    #define INDEX_REF *
    #define FROM_FORTRAN
    #define TO_FORTRAN
#endif
  using data_view_t = gridtools::data_view_tt<Value, Rank>;

public:

    IndexView(data_view_t data_view) : gt_data_view_(data_view) {
      size_ = gt_data_view_.storage_info().size();
    }

    template < typename... Coords >
    Index
    ATLAS_HOST_DEVICE
    operator()(Coords... c) {
        assert(sizeof...(Coords) == Rank);
        return INDEX_REF( &gt_data_view_(c...) );
    }

    template <typename... Coords, typename = typename boost::enable_if_c<(sizeof...(Coords) == Rank), int>::type>
    ATLAS_HOST_DEVICE
    Value const operator()(Coords... c) const {
      return gt_data_view_(c...) FROM_FORTRAN;
    }

    size_t size() const { return size_; }

    void dump(std::ostream& os) const;

private:
    data_view_t gt_data_view_;
    size_t size_;

#undef INDEX_REF
#undef TO_FORTRAN
#undef FROM_FORTRAN
};

//------------------------------------------------------------------------------------------------------

} // namespace array
} // namespace atlas
