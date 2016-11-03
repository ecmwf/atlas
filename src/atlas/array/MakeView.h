#ifndef atlas_MakeArray_h
#define atlas_MakeArray_h

#include "atlas/array/Array.h"
#include "atlas/array/ArrayView.h"
#include "atlas/array/IndexView.h"

#include "atlas/array/GridToolsTraits.h"
//------------------------------------------------------------------------------

namespace atlas {
namespace array {

#ifdef ATLAS_HAVE_GRIDTOOLS_STORAGE

template <typename Value, unsigned int NDims, bool ReadOnly = false>
static data_view_tt<Value, NDims>
make_gt_host_view(const Array& data) {
  typedef gridtools::storage_traits<BACKEND>::storage_info_t<0, NDims> storage_info_ty;
  typedef gridtools::storage_traits<BACKEND>::data_store_t<Value, storage_info_ty> data_store_t;

//      Array::data_view_t<Value, NDims, ReadOnly>
  data_store_t* ds = reinterpret_cast<data_store_t*>(((ArrayT<Value>&)data).data());

  return gridtools::make_host_view(*ds);
}


template <typename Value, unsigned int NDims, bool ReadOnly = false>
static ArrayView<Value, NDims>
make_host_view(const Array& data) {
  typedef gridtools::storage_traits<BACKEND>::storage_info_t<0, NDims> storage_info_ty;
  typedef gridtools::storage_traits<BACKEND>::data_store_t<Value, storage_info_ty> data_store_t;

//      Array::data_view_t<Value, NDims, ReadOnly>
  data_store_t* ds = reinterpret_cast<data_store_t*>(((ArrayT<Value>&)data).data());

  return ArrayView<Value, NDims>(gridtools::make_host_view(*ds));
}

#ifdef __CUDACC__
template <typename Value, unsigned int NDims, bool ReadOnly = false>
static data_view_tt<Value, NDims>
make_gt_device_view(const Array& data) {
  typedef gridtools::storage_traits<BACKEND>::storage_info_t<0, NDims> storage_info_ty;
  typedef gridtools::storage_traits<BACKEND>::data_store_t<Value, storage_info_ty> data_store_t;

//      Array::data_view_t<Value, NDims, ReadOnly>
  data_store_t* ds = reinterpret_cast<data_store_t*>(((ArrayT<Value>&)data).data());

  return gridtools::make_device_view(*ds);
}


template <typename Value, unsigned int NDims, bool ReadOnly = false>
static ArrayView<Value, NDims>
make_device_view(const Array& data) {
  typedef gridtools::storage_traits<BACKEND>::storage_info_t<0, NDims> storage_info_ty;
  typedef gridtools::storage_traits<BACKEND>::data_store_t<Value, storage_info_ty> data_store_t;

//      Array::data_view_t<Value, NDims, ReadOnly>
  data_store_t* ds = reinterpret_cast<data_store_t*>(((ArrayT<Value>&)data).data());

  return ArrayView<Value, NDims>(gridtools::make_device_view(*ds));
}
#endif

template <typename Value, unsigned int NDims, bool ReadOnly = false>
static IndexView<Value, NDims>
make_host_indexview(const Array& array) {
  typedef gridtools::storage_traits<BACKEND>::storage_info_t<0, NDims> storage_info_ty;
  typedef gridtools::storage_traits<BACKEND>::data_store_t<Value, storage_info_ty> data_store_t;

  // array::DataType<Value>::create().kind() == data.kind()
  // data.rank()

  data_store_t* ds = reinterpret_cast<data_store_t*>(((ArrayT<Value>&)array).data());

  return IndexView<Value, NDims>(gridtools::make_host_view(*ds));
}


// Old implementation
#else

template <typename Value, unsigned int NDims, bool ReadOnly = false>
static ArrayView<Value, NDims>
make_host_view(const Array& array) {
  return ArrayView<Value, NDims>((const Value*)(array.data()),array.shape());
}

template <typename Value, unsigned int NDims, bool ReadOnly = false>
static IndexView<Value, NDims>
make_host_indexview(const Array& array) {
  return IndexView<Value,NDims>( (Value*)(array.data()),array.shape().data() );
}

#endif

template <typename Value, unsigned int NDims, bool ReadOnly = false>
static ArrayView<Value, NDims>
make_view(const Array& array) {
    return make_host_view<Value, NDims, ReadOnly>(array);
}

// --------------------------------------------------------------------------------------------



static void make_host_indexview(const Array& array) {

}

template <typename Value, unsigned int NDims, bool ReadOnly = false>
static IndexView<Value, NDims>
make_indexview(const Array& array) {
  return make_host_indexview<Value,NDims>(array);
}

// template <typename Value, int NDims, bool ReadOnly = false>
// static IndexView<Value, NDims>
// make_indexview(const Value* data, const size_t idx[], const size_t shape[]) {
//   NOTIMP;
//   return IndexView<Value,NDims>();
// }

// template <>
// static IndexView<idx_t,1>
// make_indexview(const idx_t* data, const size_t idx[], const size_t shape[]) {
//   return IndexView<idx_t,1>(const_cast<idx_t*>(data+idx[0]),shape);
// }

}
}

#endif
