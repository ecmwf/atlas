/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#pragma once

#include <cassert>
#include <cstddef>
#include <cstring>

#include "atlas/array/Array.h"
#include "atlas/array/ArrayUtil.h"
#include "atlas/array/ArrayViewDefs.h"
#include "atlas/array/LocalView.h"
#include "atlas/array/gridtools/GridToolsTraits.h"
#include "atlas/library/config.h"

//------------------------------------------------------------------------------------------------------

namespace atlas {
namespace array {


template <typename Value, int Rank, Intent AccessMode = Intent::ReadWrite>
class ArrayView {
public:
    // -- Type definitions
    using value_type = typename remove_const<Value>::type;
    using return_type =
        typename std::conditional<( AccessMode == Intent::ReadWrite ), value_type, value_type const>::type;

    static constexpr Intent ACCESS{AccessMode};
    static constexpr int RANK{Rank};

    using data_view_t = gridtools::data_view_tt<value_type, Rank, gridtools::get_access_mode( AccessMode )>;

private:
    using slicer_t       = typename helpers::ArraySlicer<ArrayView<Value, Rank, AccessMode>>;
    using const_slicer_t = typename helpers::ArraySlicer<const ArrayView<Value, Rank, AccessMode>>;

    template <typename... Args>
    struct slice_t {
        using type = typename slicer_t::template Slice<Args...>::type;
    };

    template <typename... Args>
    struct const_slice_t {
        using type = typename const_slicer_t::template Slice<Args...>::type;
    };

public:
    ATLAS_HOST_DEVICE
    ArrayView( const ArrayView& other );
    ArrayView( data_view_t data_view, const Array& array );

    value_type* data() { return gt_data_view_.data(); }
    value_type const* data() const { return gt_data_view_.data(); }

    template <typename... Coords, typename = typename boost::enable_if_c<( sizeof...( Coords ) == Rank ), int>::type>
    ATLAS_HOST_DEVICE return_type& operator()( Coords... c ) {
        assert( sizeof...( Coords ) == Rank );
        return gt_data_view_( c... );
    }

    template <typename... Coords, typename = typename boost::enable_if_c<( sizeof...( Coords ) == Rank ), int>::type>
    ATLAS_HOST_DEVICE value_type const& operator()( Coords... c ) const {
        assert( sizeof...( Coords ) == Rank );
        return gt_data_view_( c... );
    }

    template <typename Int, bool EnableBool = true>
    ATLAS_HOST_DEVICE typename std::enable_if<( Rank == 1 && EnableBool ), const value_type&>::type operator[](
        Int idx ) const {
        return gt_data_view_( idx );
    }

    template <typename Int, bool EnableBool = true>
    ATLAS_HOST_DEVICE typename std::enable_if<( Rank == 1 && EnableBool ), value_type&>::type operator[]( Int idx ) {
        check_bounds( idx );
        return gt_data_view_( idx );
    }

    template <unsigned int Dim>
    ATLAS_HOST_DEVICE idx_t shape() const {
        return gt_data_view_.template length<Dim>();
    }

    ATLAS_HOST_DEVICE
    data_view_t& data_view() { return gt_data_view_; }
    ATLAS_HOST_DEVICE
    data_view_t const& data_view() const { return gt_data_view_; }

    template <unsigned int Dim>
    ATLAS_HOST_DEVICE idx_t stride() const {
        return gt_data_view_.storage_info().template stride<Dim>();
    }

    static constexpr idx_t rank() { return Rank; }
    idx_t size() const { return size_; }
    bool valid() const;

    bool contiguous() const { return ( size_ == shape_[0] * strides_[0] ? true : false ); }

    void dump( std::ostream& os ) const;

    void assign( const value_type& value );

    void assign( const std::initializer_list<value_type>& );

    const idx_t* strides() const { return strides_; }

    const idx_t* shape() const { return shape_; }

    template <typename Int>
    idx_t shape( Int idx ) const {
        return shape_[idx];
    }

    template <typename Int>
    idx_t stride( Int idx ) const {
        return strides_[idx];
    }

    template <typename... Args>
    typename slice_t<Args...>::type slice( Args... args ) {
        return slicer_t( *this ).apply( args... );
    }

    template <typename... Args>
    typename const_slice_t<Args...>::type slice( Args... args ) const {
        return const_slicer_t( *this ).apply( args... );
    }

private:
    data_view_t gt_data_view_;
    idx_t shape_[Rank];
    idx_t strides_[Rank];
    idx_t size_;
    ArrayDataStore const* data_store_orig_;
    Array const* array_;
};

//------------------------------------------------------------------------------------------------------

}  // namespace array
}  // namespace atlas
