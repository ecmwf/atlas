#ifndef atlas_MakeArray_h
#define atlas_MakeArray_h

#include "atlas/array/Array.h"
#include "atlas/array/ArrayView.h"

#include "atlas/array/GridToolsTraits.h"
//------------------------------------------------------------------------------

namespace atlas {
namespace array {

template <typename Value, unsigned int NDims, bool ReadOnly = false>
static data_view_tt<Value, NDims>
make_gt_host_view(Array* data) {
  typedef gridtools::storage_traits<BACKEND>::storage_info_t<0, NDims> storage_info_ty;
  typedef gridtools::storage_traits<BACKEND>::data_store_t<Value, storage_info_ty> data_store_t;

//      Array::data_view_t<Value, NDims, ReadOnly>
  data_store_t* ds = reinterpret_cast<data_store_t*>(((ArrayT<Value>*)data)->data());

  return gridtools::make_host_view(*ds);
}


template <typename Value, unsigned int NDims, bool ReadOnly = false>
static ArrayView<Value, NDims>
make_host_view(Array* data) {
  typedef gridtools::storage_traits<BACKEND>::storage_info_t<0, NDims> storage_info_ty;
  typedef gridtools::storage_traits<BACKEND>::data_store_t<Value, storage_info_ty> data_store_t;

//      Array::data_view_t<Value, NDims, ReadOnly>
  data_store_t* ds = reinterpret_cast<data_store_t*>(((ArrayT<Value>*)data)->data());

  return ArrayView<Value, NDims>(gridtools::make_host_view(*ds));
}

template <typename Value, unsigned int NDims, bool ReadOnly = false>
static data_view_tt<Value, NDims>
make_gt_device_view(Array* data) {
  typedef gridtools::storage_traits<BACKEND>::storage_info_t<0, NDims> storage_info_ty;
  typedef gridtools::storage_traits<BACKEND>::data_store_t<Value, storage_info_ty> data_store_t;

//      Array::data_view_t<Value, NDims, ReadOnly>
  data_store_t* ds = reinterpret_cast<data_store_t*>(((ArrayT<Value>*)data)->data());

  return gridtools::make_device_view(*ds);
}


template <typename Value, unsigned int NDims, bool ReadOnly = false>
static ArrayView<Value, NDims>
make_device_view(Array* data) {
  typedef gridtools::storage_traits<BACKEND>::storage_info_t<0, NDims> storage_info_ty;
  typedef gridtools::storage_traits<BACKEND>::data_store_t<Value, storage_info_ty> data_store_t;

//      Array::data_view_t<Value, NDims, ReadOnly>
  data_store_t* ds = reinterpret_cast<data_store_t*>(((ArrayT<Value>*)data)->data());

  return ArrayView<Value, NDims>(gridtools::make_device_view(*ds));
}


template <typename Value, unsigned int NDims, bool ReadOnly = false>
static ArrayView<Value, NDims>
make_view(Array* data) {
    return make_host_view<Value, NDims, ReadOnly>(data);
}


}
}

#endif
