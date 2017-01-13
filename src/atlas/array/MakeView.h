#pragma once

#include "atlas/array/array_fwd.h"

#include "atlas/array/Array.h"
#include "atlas/array/ArrayView.h"
#include "atlas/array/IndexView.h"
#include "atlas/array/StorageView.h"

#include "atlas/array/gridtools/GridToolsTraits.h"
//------------------------------------------------------------------------------

namespace atlas {
namespace array {

namespace impl {
    template<typename Value, unsigned NDims>
    inline static void check_metadata(const Array& array)
    {
        if(array.rank() != NDims ) {
            std::stringstream err;
            err << "Number of dimensions do not match: template argument " << NDims << " expected to be " << array.rank();
            throw eckit::BadParameter(err.str(), Here());
        }
        if(array.datatype() != array::DataType::create<Value>() ) {
            std::stringstream err;
            err << "Data Type does not match: template argument expected to be " << array.datatype().str();
            throw eckit::BadParameter(err.str(), Here());
        }
    }
}

#ifdef ATLAS_HAVE_GRIDTOOLS_STORAGE

template <typename Value, unsigned int NDims, bool ReadOnly>
inline static data_view_tt<Value, NDims>
make_gt_host_view(const Array& array) {

  typedef gridtools::storage_traits<BACKEND>::storage_info_t<0, NDims> storage_info_ty;
  typedef gridtools::storage_traits<BACKEND>::data_store_t<Value, storage_info_ty> data_store_t;

  data_store_t* ds = reinterpret_cast<data_store_t*>(const_cast<void*>(array.storage()));

  return gridtools::make_host_view(*ds);
}


template <typename Value, unsigned int NDims, bool ReadOnly>
inline ArrayView<Value, NDims>
make_host_view(const Array& array) {
  impl::check_metadata<Value, NDims>(array);
  return ArrayView<Value, NDims>(make_gt_host_view<Value, NDims>(array), array);
}

template <typename Value>
inline StorageView<Value>
make_host_storageview(const Array& array) {
  typedef gridtools::storage_traits<BACKEND>::storage_info_t<0, 1> storage_info_ty;
  typedef gridtools::storage_traits<BACKEND>::data_store_t<Value, storage_info_ty> data_store_t;
  data_store_t* ds = reinterpret_cast<data_store_t*>(const_cast<void*>(array.storage()));
  return StorageView<Value>(gridtools::make_host_view(*ds),array.size(),array.contiguous());
}

#ifdef ENABLE_GPU
template <typename Value, unsigned int NDims, bool ReadOnly = false>
inline static data_view_tt<Value, NDims>
make_gt_device_view(const Array& array) {
  typedef gridtools::storage_traits<BACKEND>::storage_info_t<0, NDims> storage_info_ty;
  typedef gridtools::storage_traits<BACKEND>::data_store_t<Value, storage_info_ty> data_store_t;

  data_store_t* ds = reinterpret_cast<data_store_t*>(const_cast<void*>(array.storage()));

  return gridtools::make_device_view(*ds);
}

template <typename Value, unsigned int NDims, bool ReadOnly>
inline ArrayView<Value, NDims>
make_device_view(const Array& array) {
  impl::check_metadata<Value, NDims>(array);
  return ArrayView<Value, NDims>(make_gt_device_view<Value, NDims>(array), array);
}


template <typename Value>
inline StorageView<Value>
make_device_storageview(const Array& array) {
  return StorageView<Value>(make_gt_device_view<Value, 1>(array),array.size(),array.contiguous());
}

#else
template <typename Value, unsigned int NDims, bool ReadOnly>
inline ArrayView<Value, NDims>
make_device_view(const Array& array) {
  impl::check_metadata<Value, NDims>(array);
  return ArrayView<Value, NDims>(make_gt_host_view<Value, NDims>(array), array);
}

template <typename Value>
inline StorageView<Value>
make_device_storageview(const Array& array) {
  return StorageView<Value>(make_gt_host_view<Value, 1>(array),array.size(),array.contiguous());
}

#endif

template <typename Value, unsigned int NDims, bool ReadOnly>
inline IndexView<Value, NDims>
make_host_indexview(const Array& array) {
  typedef gridtools::storage_traits<BACKEND>::storage_info_t<0, NDims> storage_info_ty;
  typedef gridtools::storage_traits<BACKEND>::data_store_t<Value, storage_info_ty> data_store_t;

  data_store_t* ds = reinterpret_cast<data_store_t*>(const_cast<void*>(array.storage()));

  return IndexView<Value, NDims>(gridtools::make_host_view(*ds));
}


// Old implementation
#else

template <typename Value, unsigned int NDims, bool ReadOnly>
inline ArrayView<Value, NDims>
make_host_view(const Array& array) {
  return ArrayView<Value, NDims>((const Value*)(array.storage()),array.shape());
}

template <typename Value, unsigned int NDims, bool ReadOnly>
inline ArrayView<Value, NDims>
make_device_view(const Array& array) {
  return make_host_view<Value,NDims,ReadOnly>(array);
}

template <typename Value, unsigned int NDims, bool ReadOnly>
inline IndexView<Value, NDims>
make_host_indexview(const Array& array) {
  return IndexView<Value,NDims>( (Value*)(array.storage()),array.shape().data() );
}

template <typename Value>
inline StorageView<Value>
make_host_storageview(const Array& array) {
  return StorageView<Value>(const_cast<Array&>(array).storage(),array.size(),array.contiguous());
}


template <typename Value>
inline StorageView<Value>
make_device_storageview(const Array& array) {
  return make_host_storageview<Value>(array);
}

#endif

// --------------------------------------------------------------------------------------------

template <typename Value, unsigned int NDims, bool ReadOnly>
inline IndexView<Value, NDims>
make_indexview(const Array& array) {
  impl::check_metadata<Value, NDims>(array);
  return make_host_indexview<Value,NDims>(array);
}

template <typename Value, unsigned int NDims, bool ReadOnly>
inline ArrayView<Value, NDims>
make_view(const Array& array) {
    impl::check_metadata<Value, NDims>(array);

    return make_host_view<Value, NDims, ReadOnly>(array);
}

template <typename Value>
inline StorageView<Value>
make_storageview(const Array& array) {
    return make_host_storageview<Value>(array);
}


}
}
