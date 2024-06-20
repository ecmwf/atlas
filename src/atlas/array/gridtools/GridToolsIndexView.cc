/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <ostream>
#include <stdexcept>

#include "atlas/array.h"
#include "atlas/array/IndexView.h"
#include "atlas/field/Field.h"
#include "atlas/runtime/Exception.h"

//------------------------------------------------------------------------------------------------------

namespace atlas {
namespace array {

namespace detail {

template <typename T, size_t Rank>
struct host_device_array { // Copied from GridToolsArrayView.cc
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
struct StoragePropBuilder; // Copied from GridToolsArrayView.cc

template <typename StorageInfo, int Rank, std::size_t... I>
struct StoragePropBuilder<StorageInfo, Rank, std::integer_sequence<std::size_t, I...>> { // Copied from GridToolsArrayView.cc
    static host_device_array<ArrayStrides::value_type, Rank> buildStrides(const StorageInfo& storage_info) {
        return host_device_array<ArrayStrides::value_type, Rank>{
            (ArrayStrides::value_type)storage_info.template stride<I>()...};
    }
    static host_device_array<ArrayShape::value_type, Rank> buildShapes(const StorageInfo& storage_info) {
        return host_device_array<ArrayShape::value_type, Rank>{storage_info.template total_length<I>()...};
    }
};

}


template <typename Value, int Rank>
IndexView<Value, Rank>::IndexView(const Array& array, bool _device_view):
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
            detail::StoragePropBuilder<storage_info_ty, Rank, std::make_integer_sequence<std::size_t, Rank>>::buildStrides(
                storage_info_);
        auto shapet =
            detail::StoragePropBuilder<storage_info_ty, Rank, std::make_integer_sequence<std::size_t, Rank>>::buildShapes(
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



    // IndexView::IndexView(data_view_t data_view): gt_data_view_(data_view) {
    //     if (gt_data_view_.valid()) {
    //         size_ = gt_data_view_.storage_info().total_length();
    //     }
    //     else {
    //         size_ = 0;
    //     }



    // if (gt_data_view_.valid()) {
    //     constexpr static unsigned int ndims = data_view_t::data_store_t::storage_info_t::ndims;

    //     using storage_info_ty = gridtools::storage_traits::storage_info_t<0, ndims>;
    //     using data_store_t    = gridtools::storage_traits::data_store_t<value_type, storage_info_ty>;

    //     auto storage_info_ =
    //         *((reinterpret_cast<data_store_t*>(const_cast<void*>(array.storage())))->get_storage_info_ptr());

    //     auto stridest =
    //         StoragePropBuilder<storage_info_ty, Rank, std::make_integer_sequence<std::size_t, Rank>>::buildStrides(
    //             storage_info_);
    //     auto shapet =
    //         StoragePropBuilder<storage_info_ty, Rank, std::make_integer_sequence<std::size_t, Rank>>::buildShapes(
    //             storage_info_);


    //     for (int i = 0; i < Rank; ++i) {
    //         strides_[i] = stridest[i];
    //         shape_[i]   = shapet[i];
    //     }

    //     size_ = storage_info_.total_length();
    // }
    // else {
    //     std::fill_n(shape_, Rank, 0);
    //     std::fill_n(strides_, Rank, 0);

    //     size_ = 0;
    // }

    // }


//------------------------------------------------------------------------------------------------------

template <typename Value, int Rank>
void IndexView<Value, Rank>::dump(std::ostream& os) const {
    ATLAS_NOTIMPLEMENTED;
}


template <typename Value, int Rank>
LocalIndexView<Value, Rank>::LocalIndexView(Value* data, const idx_t shape[1]): data_(const_cast<Value*>(data)) {
    strides_[0] = 1;
    shape_[0]   = shape[0];
}

template <typename Value, int Rank>
LocalIndexView<Value, Rank>::LocalIndexView(Value* data, const idx_t shape[1], const idx_t strides[1]):
    data_(const_cast<Value*>(data)) {
    strides_[0] = strides[0];
    shape_[0]   = shape[0];
}


//------------------------------------------------------------------------------------------------------

}  // namespace array
}  // namespace atlas

//------------------------------------------------------------------------------------------------------
// Explicit instantiation
namespace atlas {
namespace array {

#define EXPLICIT_TEMPLATE_INSTANTIATION_TYPE_RANK(TYPE, RANK) \
    template class IndexView<TYPE, RANK>;                     \
    template class IndexView<const TYPE, RANK>;               \
    template class LocalIndexView<TYPE, RANK>;                \
    template class LocalIndexView<const TYPE, RANK>;

#define EXPLICIT_TEMPLATE_INSTANTIATION(RANK)            \
    EXPLICIT_TEMPLATE_INSTANTIATION_TYPE_RANK(int, RANK) \
    EXPLICIT_TEMPLATE_INSTANTIATION_TYPE_RANK(long, RANK)

EXPLICIT_TEMPLATE_INSTANTIATION(1)
EXPLICIT_TEMPLATE_INSTANTIATION(2)

#undef EXPLICIT_TEMPLATE_INSTANTIATION_TYPE_RANK
#undef EXPLICIT_TEMPLATE_INSTANTIATION

}  // namespace array
}  // namespace atlas
