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
#include "atlas/array/helpers/ArrayCopier.h"
#include "atlas/array/helpers/ArrayInitializer.h"
#include "atlas/array/helpers/ArrayWriter.h"
#include "atlas/runtime/Exception.h"

#define ENABLE_IF_NON_CONST \
    template <bool EnableBool, typename std::enable_if<(!std::is_const<Value>::value && EnableBool), int>::type*>

#define ENABLE_IF_CONST \
    template <bool EnableBool, typename std::enable_if<(std::is_const<Value>::value && EnableBool), int>::type*>

namespace atlas {
namespace array {


template <typename T, size_t Rank>
struct host_device_array {
    ATLAS_HOST_DEVICE host_device_array(std::initializer_list<T> list) {
        size_t i(0);
        for (const T v : list) {
            data_[i++] = v;
        }
    }
    ATLAS_HOST_DEVICE ~host_device_array() {}

    ATLAS_HOST_DEVICE const T* data() const { return data_; }

    T operator[](int i) const { return data_[i]; }

    T data_[Rank];
};


template <typename StorageInfo, int Rank, typename IS>
struct StoragePropBuilder;

template <typename StorageInfo, int Rank, std::size_t... I>
struct StoragePropBuilder<StorageInfo, Rank, std::integer_sequence<std::size_t, I...>> {
    static host_device_array<ArrayStrides::value_type, Rank> buildStrides(const StorageInfo& storage_info) {
        return host_device_array<ArrayStrides::value_type, Rank>{
            (ArrayStrides::value_type)storage_info.template stride<I>()...};
    }
    static host_device_array<ArrayShape::value_type, Rank> buildShapes(const StorageInfo& storage_info) {
        return host_device_array<ArrayShape::value_type, Rank>{storage_info.template total_length<I>()...};
    }
};

template <typename Value, int Rank>
ArrayView<Value, Rank>::ArrayView(const ArrayView& other):
    gt_data_view_(other.gt_data_view_),
    data_store_orig_(other.data_store_orig_),
    array_(other.array_),
    is_device_view_(other.is_device_view_) {
    std::memcpy(shape_, other.shape_, sizeof(ArrayShape::value_type) * Rank);
    std::memcpy(strides_, other.strides_, sizeof(ArrayStrides::value_type) * Rank);
    size_ = other.size_;
    // TODO: check compatibility
}

template <typename Value, int Rank>
ArrayView<Value, Rank>::ArrayView(const Array& array, bool _device_view):
    gt_data_view_(_device_view ? gridtools::make_gt_device_view<Value, Rank>(array)
                               : gridtools::make_gt_host_view<Value, Rank>(array)),
    data_store_orig_(&array.data_store()),
    array_(&array),
    is_device_view_(_device_view) {
    if (gt_data_view_.valid()) {
        constexpr static unsigned int ndims = data_view_t::data_store_t::storage_info_t::ndims;

        using storage_info_ty = gridtools::storage_traits::storage_info_t<0, ndims>;
        using data_store_t    = gridtools::storage_traits::data_store_t<value_type, storage_info_ty>;

        auto storage_info_ =
            *((reinterpret_cast<data_store_t*>(const_cast<void*>(array.storage())))->get_storage_info_ptr());

        auto stridest =
            StoragePropBuilder<storage_info_ty, Rank, std::make_integer_sequence<std::size_t, Rank>>::buildStrides(
                storage_info_);
        auto shapet =
            StoragePropBuilder<storage_info_ty, Rank, std::make_integer_sequence<std::size_t, Rank>>::buildShapes(
                storage_info_);


        for (int i = 0; i < Rank; ++i) {
            strides_[i] = stridest[i];
            shape_[i]   = shapet[i];
        }

        size_ = storage_info_.total_length();
    }
    else {
        std::fill_n(shape_, Rank, 0);
        std::fill_n(strides_, Rank, 0);

        size_ = 0;
    }
}

template <typename Value, int Rank>
bool ArrayView<Value, Rank>::valid() const {
    return gt_data_view_.valid() && (&array_->data_store() == data_store_orig_);
}


template <typename Value, int Rank>
ENABLE_IF_NON_CONST void ArrayView<Value, Rank>::assign(const value_type& value) {
    helpers::array_assigner<Value, Rank>::apply(*this, value);
}

//------------------------------------------------------------------------------------------------------

template <typename Value, int Rank>
ENABLE_IF_NON_CONST void ArrayView<Value, Rank>::assign(const std::initializer_list<value_type>& list) {
    helpers::array_assigner<Value, Rank>::apply(*this, list);
}

//------------------------------------------------------------------------------------------------------


template <typename Value, int Rank>
ENABLE_IF_NON_CONST void ArrayView<Value, Rank>::assign(const ArrayView& other) {
    helpers::array_copier<Value, Rank>::apply(other, *this);
}

//------------------------------------------------------------------------------------------------------

template <typename Value, int Rank>
void ArrayView<Value, Rank>::dump(std::ostream& os) const {
    os << "size: " << size() << " , values: ";
    os << "[ ";
    helpers::array_writer::apply(*this, os);
    os << " ]";
}

//------------------------------------------------------------------------------------------------------

}  // namespace array
}  // namespace atlas

//-----------------------------------------------------------------------
// Explicit instantiation
namespace atlas {
namespace array {
#define EXPLICIT_TEMPLATE_INSTANTIATION(Rank)                                                                         \
    template class ArrayView<int, Rank>;                                                                              \
    template class ArrayView<int const, Rank>;                                                                        \
    template class ArrayView<long, Rank>;                                                                             \
    template class ArrayView<long const, Rank>;                                                                       \
    template class ArrayView<long unsigned, Rank>;                                                                    \
    template class ArrayView<long unsigned const, Rank>;                                                              \
    template class ArrayView<float, Rank>;                                                                            \
    template class ArrayView<float const, Rank>;                                                                      \
    template class ArrayView<double, Rank>;                                                                           \
    template class ArrayView<double const, Rank>;                                                                     \
                                                                                                                      \
    template void ArrayView<int, Rank>::assign<true, nullptr>(int const&);                                            \
    template void ArrayView<long, Rank>::assign<true, nullptr>(long const&);                                          \
    template void ArrayView<float, Rank>::assign<true, nullptr>(float const&);                                        \
    template void ArrayView<double, Rank>::assign<true, nullptr>(double const&);                                      \
    template void ArrayView<long unsigned, Rank>::assign<true, nullptr>(long unsigned const&);                        \
    template void ArrayView<int, Rank>::assign<true, nullptr>(std::initializer_list<int> const&);                     \
    template void ArrayView<long, Rank>::assign<true, nullptr>(std::initializer_list<long> const&);                   \
    template void ArrayView<float, Rank>::assign<true, nullptr>(std::initializer_list<float> const&);                 \
    template void ArrayView<double, Rank>::assign<true, nullptr>(std::initializer_list<double> const&);               \
    template void ArrayView<long unsigned, Rank>::assign<true, nullptr>(std::initializer_list<long unsigned> const&); \
    template void ArrayView<int, Rank>::assign<true, nullptr>(ArrayView<int, Rank> const&);                           \
    template void ArrayView<long, Rank>::assign<true, nullptr>(ArrayView<long, Rank> const&);                         \
    template void ArrayView<float, Rank>::assign<true, nullptr>(ArrayView<float, Rank> const&);                       \
    template void ArrayView<double, Rank>::assign<true, nullptr>(ArrayView<double, Rank> const&);


// For each Rank in [1..9]
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
}  // namespace array
}  // namespace atlas
