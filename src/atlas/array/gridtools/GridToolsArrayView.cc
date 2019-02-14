/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/array/gridtools/GridToolsArrayView.h"
#include <cassert>
#include "atlas/array/gridtools/GridToolsArrayHelpers.h"
#include "atlas/array/helpers/ArrayAssigner.h"
#include "atlas/array/helpers/ArrayInitializer.h"
#include "atlas/array/helpers/ArrayWriter.h"
#include "atlas/runtime/Exception.h"

namespace atlas {
namespace array {

template <typename T, size_t Rank>
struct host_device_array {
    ATLAS_HOST_DEVICE host_device_array( std::initializer_list<T> list ) {
        size_t i( 0 );
        for ( const T v : list ) {
            data_[i++] = v;
        }
    }
    ATLAS_HOST_DEVICE ~host_device_array() {}

    ATLAS_HOST_DEVICE const T* data() const { return data_; }

    T operator[]( int i ) const { return data_[i]; }

    T data_[Rank];
};

template <typename Value, int Rank, Intent AccessMode>
ArrayView<Value, Rank, AccessMode>::ArrayView( const ArrayView& other ) :
    gt_data_view_( other.gt_data_view_ ),
    data_store_orig_( other.data_store_orig_ ),
    array_( other.array_ ) {
    std::memcpy( shape_, other.shape_, sizeof( ArrayShape::value_type ) * Rank );
    std::memcpy( strides_, other.strides_, sizeof( ArrayStrides::value_type ) * Rank );
    size_ = other.size_;
    // TODO: check compatibility
}

template <typename Value, int Rank, Intent AccessMode>
ArrayView<Value, Rank, AccessMode>::ArrayView( data_view_t data_view, const Array& array ) :
    gt_data_view_( data_view ),
    data_store_orig_( array.data_store() ),
    array_( &array ) {
    if ( data_view.valid() ) {
        using seq = ::gridtools::apply_gt_integer_sequence<::gridtools::make_gt_integer_sequence<int, Rank>>;

        constexpr static unsigned int ndims = data_view_t::data_store_t::storage_info_t::ndims;

        using storage_info_ty = gridtools::storage_traits::storage_info_t<0, ndims>;
        using data_store_t    = gridtools::storage_traits::data_store_t<value_type, storage_info_ty>;

        auto storage_info_ =
            *( ( reinterpret_cast<data_store_t*>( const_cast<void*>( array.storage() ) ) )->get_storage_info_ptr() );

        auto stridest = seq::template apply<
            host_device_array<ArrayStrides::value_type, Rank>,
            atlas::array::gridtools::get_stride_component<ArrayStrides::value_type>::template get_component>(
            &( storage_info_ ) );
        auto shapet = seq::template apply<
            host_device_array<ArrayShape::value_type, Rank>,
            atlas::array::gridtools::get_shape_component<ArrayStrides::value_type>::template get_component>(
            &( storage_info_ ) );

        for ( int i = 0; i < Rank; ++i ) {
            strides_[i] = stridest[i];
            shape_[i]   = shapet[i];
        }

        size_ = storage_info_.total_length();
    }
    else {
        std::fill_n( shape_, Rank, 0 );
        std::fill_n( strides_, Rank, 0 );

        size_ = 0;
    }
}

template <typename Value, int Rank, Intent AccessMode>
bool ArrayView<Value, Rank, AccessMode>::valid() const {
    return gt_data_view_.valid() && ( array_->data_store() == data_store_orig_ );
}

template <typename Value, int Rank, Intent AccessMode>
void ArrayView<Value, Rank, AccessMode>::assign( const value_type& value ) {
    helpers::array_assigner<Value, Rank>::apply( *this, value );
}

template <typename Value, int Rank, Intent AccessMode>
void ArrayView<Value, Rank, AccessMode>::assign( const std::initializer_list<value_type>& list ) {
    ATLAS_ASSERT( list.size() == size_ );
    helpers::array_assigner<Value, Rank>::apply( *this, list );
}

template <typename Value, int Rank, Intent AccessMode>
void ArrayView<Value, Rank, AccessMode>::dump( std::ostream& os ) const {
    os << "size: " << size() << " , values: ";
    os << "[ ";
    helpers::array_writer::apply( *this, os );
    os << " ]";
}

//------------------------------------------------------------------------------------------------------

}  // namespace array
}  // namespace atlas

//-----------------------------------------------------------------------
// Explicit instantiation
namespace atlas {
namespace array {
#define EXPLICIT_TEMPLATE_INSTANTIATION( Rank )                       \
    template class ArrayView<int, Rank, Intent::ReadOnly>;            \
    template class ArrayView<int, Rank, Intent::ReadWrite>;           \
    template class ArrayView<long, Rank, Intent::ReadOnly>;           \
    template class ArrayView<long, Rank, Intent::ReadWrite>;          \
    template class ArrayView<long unsigned, Rank, Intent::ReadOnly>;  \
    template class ArrayView<long unsigned, Rank, Intent::ReadWrite>; \
    template class ArrayView<float, Rank, Intent::ReadOnly>;          \
    template class ArrayView<float, Rank, Intent::ReadWrite>;         \
    template class ArrayView<double, Rank, Intent::ReadOnly>;         \
    template class ArrayView<double, Rank, Intent::ReadWrite>;

// For each NDims in [1..9]
EXPLICIT_TEMPLATE_INSTANTIATION( 1 )
EXPLICIT_TEMPLATE_INSTANTIATION( 2 )
EXPLICIT_TEMPLATE_INSTANTIATION( 3 )
EXPLICIT_TEMPLATE_INSTANTIATION( 4 )
EXPLICIT_TEMPLATE_INSTANTIATION( 5 )
EXPLICIT_TEMPLATE_INSTANTIATION( 6 )
EXPLICIT_TEMPLATE_INSTANTIATION( 7 )
EXPLICIT_TEMPLATE_INSTANTIATION( 8 )
EXPLICIT_TEMPLATE_INSTANTIATION( 9 )

#undef EXPLICIT_TEMPLATE_INSTANTIATION
}  // namespace array
}  // namespace atlas
