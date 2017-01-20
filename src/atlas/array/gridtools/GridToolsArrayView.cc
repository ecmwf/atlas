/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "atlas/array/gridtools/GridToolsArrayView.h"
#include "atlas/array/gridtools/GridToolsArrayHelpers.h"
#include "eckit/exception/Exceptions.h"

namespace atlas {
namespace array {

template< typename Value, int Rank >
ArrayView<Value,Rank>::ArrayView( const ArrayView& other ) :
    gt_data_view_(other.gt_data_view_) {
    std::memcpy(shape_,other.shape_,sizeof(size_t)*Rank);
    std::memcpy(strides_,other.strides_,sizeof(size_t)*Rank);
    size_ = other.size_;
    // TODO: check compatibility
}

template< typename Value, int Rank >
ArrayView<Value,Rank>::ArrayView(data_view_t data_view, const Array& array) :
    gt_data_view_(data_view) {
    using seq = ::gridtools::apply_gt_integer_sequence<typename ::gridtools::make_gt_integer_sequence<int, Rank>::type>;

    constexpr static unsigned int ndims = data_view_t::data_store_t::storage_info_t::ndims;
    data_view_t gt_host_view_ = atlas::array::gridtools::make_gt_host_view<value_type, ndims, true> ( array );

    auto stridest = seq::template apply<
        std::vector<size_t>,
        atlas::array::gridtools::get_stride_component<unsigned long, ::gridtools::static_uint<Rank> >::template get_component>(
        &(gt_host_view_.storage_info()));
    auto shapet = seq::template apply<std::vector<size_t>, atlas::array::gridtools::get_shape_component>(&(gt_host_view_.storage_info()));

    std::memcpy(strides_, &(stridest[0]), sizeof(size_t)*Rank);
    std::memcpy(shape_, &(shapet[0]), sizeof(size_t)*Rank);

    size_ = gt_host_view_.storage_info().size();
}

template< typename Value, int Rank >
void ArrayView<Value,Rank>::assign(const value_type& value) {
    ASSERT( contiguous() );
    value_type* raw_data = data();
    for( size_t j=0; j<size_; ++j ) {
        raw_data[j] = value;
    }
}

template< typename Value, int Rank >
void ArrayView<Value,Rank>::dump(std::ostream& os) const {
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
template class ArrayView<int,Rank>;\
template class ArrayView<long,Rank>;\
template class ArrayView<long unsigned,Rank>;\
template class ArrayView<float,Rank>;\
template class ArrayView<double,Rank>;\

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
