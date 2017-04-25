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

#include <cstddef>
#include <cstring>
#include <cassert>
#include "atlas/library/config.h"
#include "atlas/array/ArrayUtil.h"
#include "atlas/array/LocalView.h"
#include "atlas/array/gridtools/GridToolsTraits.h"
#include "atlas/array/gridtools/GridToolsMakeView.h"

//------------------------------------------------------------------------------------------------------

namespace atlas {
namespace array {

template< typename Value, int Rank > class ArrayView {
public:
// -- Type definitions
    using value_type = typename remove_const<Value>::type;
    using Slice = typename std::conditional<(Rank==1), value_type&, LocalView<value_type,Rank-1> >::type;
    using data_view_t = atlas::array::gridtools::data_view_tt<value_type, Rank>;

public:

    ArrayView( const ArrayView& other );
    ArrayView(data_view_t data_view, const Array& array);
    value_type* data() { return gt_data_view_.data(); }
    value_type const* data() const { return gt_data_view_.data(); }

    template < typename... Coords, typename = typename boost::enable_if_c<(sizeof...(Coords) == Rank), int>::type >
    value_type&
    ATLAS_HOST_DEVICE
    operator()(Coords... c) {
        assert(sizeof...(Coords) == Rank);
        return gt_data_view_(c...);
    }

    template <typename... Coords, typename = typename boost::enable_if_c<(sizeof...(Coords) == Rank), int>::type>
    ATLAS_HOST_DEVICE
    value_type const& operator()(Coords... c) const {
        return gt_data_view_(c...);
    }

    size_t shape(size_t idx) const { return shape_[idx]; }

    const Slice operator[](size_t i) const {
        return Slicer<Slice, Rank==1>(*this).apply(i);
    }

    Slice operator[](size_t i) {
        return Slicer<Slice, Rank==1>(*this).apply(i);
    }

    const Slice at(size_t i) const {
        if( i>= shape(0) ) throw eckit::OutOfRange(i,shape(0),Here());
        return Slicer<Slice, Rank==1>(*this).apply(i);
    }

    Slice at(size_t i) {
        if( i>= shape(0) ) throw eckit::OutOfRange(i,shape(0),Here());
        return Slicer<Slice, Rank==1>(*this).apply(i);
    }


    data_view_t& data_view() { return gt_data_view_;}

    size_t rank() const { return Rank; }
    size_t size() const { return size_; }
    bool valid() const;

    bool contiguous() const {
        return (size_ == shape_[0]*strides_[0] ? true : false);
    }

    void dump(std::ostream& os) const;

    void assign(const value_type& value);

private:

    template <typename ReturnType = Slice, bool ToScalar = false>
    struct Slicer {
        Slicer(ArrayView<value_type, Rank> const& av) : av_(av) {}
        ArrayView<value_type, Rank> const& av_;
        ReturnType apply(const size_t i) const {
          return LocalView<value_type,Rank-1>(
                  av_.data()+av_.strides_[0]*i,
                  av_.shape_+1,
                  av_.strides_+1 );
        }
    };

    template <typename ReturnType>
    struct Slicer<ReturnType, true> {
        Slicer(ArrayView<value_type, Rank> const& av) : av_(av) {}
        ArrayView<value_type, Rank> const& av_;
        ReturnType apply(const size_t i) const {
            return *(const_cast<value_type*>(av_.data()) + av_.strides_[0] * i);
        }
    };


    data_view_t gt_data_view_;
    size_t shape_[Rank];
    size_t strides_[Rank];
    size_t size_;
    ArrayDataStore const* data_store_orig_;
    Array const* array_;
};

//------------------------------------------------------------------------------------------------------

} // namespace array
} // namespace atlas
