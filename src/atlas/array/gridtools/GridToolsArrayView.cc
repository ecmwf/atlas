/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <cassert>
#include "atlas/array/gridtools/GridToolsArrayView.h"
#include "atlas/array/gridtools/GridToolsArrayHelpers.h"
#include "atlas/array/helpers/ArrayInitializer.h"
#include "atlas/array/helpers/ArrayAssigner.h"
#include "eckit/exception/Exceptions.h"

namespace atlas {
namespace array {

template< typename T >
struct private_vector {

    ATLAS_HOST_DEVICE private_vector( std::initializer_list<T> list ) {
        data_ = new T( list.size() );
        size_t i(0);
        for( const T v : list ) {
            data_[i++] = v;
        }
    }
    ATLAS_HOST_DEVICE ~private_vector() {
        delete data_;
    }

    ATLAS_HOST_DEVICE const T* data() const {
        return data_;
    }

    T* data_;
};

template< typename Value, int Rank, Intent AccessMode >
ArrayView<Value,Rank,AccessMode>::ArrayView( const ArrayView& other ) :
    gt_data_view_(other.gt_data_view_), data_store_orig_(other.data_store_orig_), array_(other.array_) {
    std::memcpy(shape_,other.shape_,sizeof(size_t)*Rank);
    std::memcpy(strides_,other.strides_,sizeof(size_t)*Rank);
    size_ = other.size_;
    // TODO: check compatibility
}

template< typename Value, int Rank, Intent AccessMode >
ArrayView<Value,Rank, AccessMode>::ArrayView(data_view_t data_view, const Array& array) :
    gt_data_view_(data_view), data_store_orig_(array.data_store()), array_(&array) {
    if(data_view.valid()) {
        using seq = ::gridtools::apply_gt_integer_sequence<typename ::gridtools::make_gt_integer_sequence<int, Rank>::type>;

        constexpr static unsigned int ndims = data_view_t::data_store_t::storage_info_t::ndims;

        using storage_info_ty = gridtools::storage_traits::storage_info_t<0, ndims>;
        using data_store_t    = gridtools::storage_traits::data_store_t<value_type, storage_info_ty>;

        auto storage_info_ = *((reinterpret_cast<data_store_t*>(const_cast<void*>(array.storage())))->get_storage_info_ptr());

        auto stridest = seq::template apply<
            private_vector<size_t>,
            atlas::array::gridtools::get_stride_component<unsigned long, ::gridtools::static_uint<Rank> >::template get_component>(
            &(storage_info_));
        auto shapet = seq::template apply<
            private_vector<size_t>,
            atlas::array::gridtools::get_shape_component>(&(storage_info_));

        std::memcpy(strides_, stridest.data(), sizeof(size_t)*Rank);
        std::memcpy(shape_, shapet.data(), sizeof(size_t)*Rank);

        size_ = storage_info_.total_length();
    }
    else {

        std::fill_n(shape_, Rank, 0 );
        std::fill_n(strides_, Rank, 0 );

        size_ = 0;
    }
}

template< typename Value, int Rank, Intent AccessMode>
bool ArrayView<Value,Rank, AccessMode>::valid() const {
    return gt_data_view_.valid() && (array_->data_store() == data_store_orig_);
}

template< typename Value, int Rank, Intent AccessMode >
void ArrayView<Value,Rank,AccessMode>::assign(const value_type& value) {
    helpers::array_assigner<Value,Rank>::apply(*this,value);
}

template <typename Value, int Rank, Intent AccessMode>
void ArrayView<Value,Rank,AccessMode>::assign(const std::initializer_list<value_type>& list) {
    ASSERT( list.size() == size_ );
    ASSERT( contiguous() );
    value_type* raw_data = data();
    size_t j(0);
    for( const value_type& v : list ) {
        raw_data[j++] = v;
    }
}

template< typename Value, int Rank, Intent AccessMode >
void ArrayView<Value,Rank,AccessMode>::dump(std::ostream& os) const {
    ASSERT( contiguous() );

    const value_type* data_ = data();
    os << "size: " << size() << " , values: ";
    os << "[ ";
    for( size_t j=0; j<size(); ++ j )
        os << data_[j] << " ";
    os << "]";
}

//------------------------------------------------------------------------------------------------------

} // namespace array
} // namespace atlas

//-----------------------------------------------------------------------
// Explicit instantiation
namespace atlas {
namespace array {
#define EXPLICIT_TEMPLATE_INSTANTIATION(Rank) \
template class ArrayView<int,Rank,Intent::ReadOnly>;\
template class ArrayView<int,Rank,Intent::ReadWrite>;\
template class ArrayView<long,Rank,Intent::ReadOnly>;\
template class ArrayView<long,Rank,Intent::ReadWrite>;\
template class ArrayView<long unsigned,Rank,Intent::ReadOnly>;\
template class ArrayView<long unsigned,Rank,Intent::ReadWrite>;\
template class ArrayView<float,Rank,Intent::ReadOnly>;\
template class ArrayView<float,Rank,Intent::ReadWrite>;\
template class ArrayView<double,Rank,Intent::ReadOnly>;\
template class ArrayView<double,Rank,Intent::ReadWrite>;\

// For each NDims in [1..9]
EXPLICIT_TEMPLATE_INSTANTIATION(1)
EXPLICIT_TEMPLATE_INSTANTIATION(2)
EXPLICIT_TEMPLATE_INSTANTIATION(3)
EXPLICIT_TEMPLATE_INSTANTIATION(4)
EXPLICIT_TEMPLATE_INSTANTIATION(5)
EXPLICIT_TEMPLATE_INSTANTIATION(6)
EXPLICIT_TEMPLATE_INSTANTIATION(7)
EXPLICIT_TEMPLATE_INSTANTIATION(8)
EXPLICIT_TEMPLATE_INSTANTIATION(9)

#undef EXPLICIT_TEMPLATE_INSTANTIATION
}
}
