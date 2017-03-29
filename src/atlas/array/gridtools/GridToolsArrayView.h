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

//#include "atlas/array/native/ArrayView_iterator.h"

//------------------------------------------------------------------------------------------------------

namespace atlas {
namespace array {

template< typename Value, int Rank >
class ArrayView
{
public:
// -- Type definitions
//  typedef ArrayView_iterator<Value,Rank>         iterator;
//  typedef ArrayView_const_iterator<Value,Rank>   const_iterator;
    typedef typename remove_const<Value>::type  value_type;

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

    LocalView<value_type,Rank-1> at(const size_t i) const {
        assert( i < shape_[0] );
        return LocalView<value_type,Rank-1>(
                const_cast<value_type*>(data())+strides_[0]*i,
                shape_+1,
                strides_+1 );
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
