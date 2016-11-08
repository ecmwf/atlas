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

namespace impl {
    template<typename Value, unsigned NDims>
    inline static void check_metadata(const Array& array)
    {
        if(array.rank() != NDims ) {
            std::stringstream err;
            err << "Number of dimensions do not match";
            throw eckit::BadParameter(err.str(), Here());
        }
        if(array.datatype() != array::DataType::create<Value>() ) {
            std::stringstream err;
            err << "Data Type does not match";
            throw eckit::BadParameter(err.str(), Here());
        }
    }
}

template <typename Value, unsigned int NDims, bool ReadOnly = false>
inline static data_view_tt<Value, NDims>
make_gt_host_view(const Array& array) {

  impl::check_metadata<Value, NDims>(array);
  typedef gridtools::storage_traits<BACKEND>::storage_info_t<0, NDims> storage_info_ty;
  typedef gridtools::storage_traits<BACKEND>::data_store_t<Value, storage_info_ty> data_store_t;

//      Array::data_view_t<Value, NDims, ReadOnly>
  data_store_t* ds = reinterpret_cast<data_store_t*>(const_cast<void*>(array.data()));

  return gridtools::make_host_view(*ds);
}


template <typename Value, unsigned int NDims, bool ReadOnly = false>
inline static ArrayView<Value, NDims>
make_host_view(const Array& array) {
  return ArrayView<Value, NDims>(make_gt_host_view<Value, NDims>(array), array.shape());
}


#ifdef __CUDACC__
template <typename Value, unsigned int NDims, bool ReadOnly = false>
inline static data_view_tt<Value, NDims>
make_gt_device_view(const Array& array) {
  impl::check_metadata<Value, NDims>(array);

  typedef gridtools::storage_traits<BACKEND>::storage_info_t<0, NDims> storage_info_ty;
  typedef gridtools::storage_traits<BACKEND>::data_store_t<Value, storage_info_ty> data_store_t;

//      Array::data_view_t<Value, NDims, ReadOnly>
  data_store_t* ds = reinterpret_cast<data_store_t*>(array.data());

  return gridtools::make_device_view(*ds);
}


template <typename Value, unsigned int NDims, bool ReadOnly = false>
inline static ArrayView<Value, NDims>
make_device_view(const Array& array) {
  return ArrayView<Value, NDims>(make_gt_host_view<Value, NDims>(array), ds->shape());
}
#endif

template <typename Value, unsigned int NDims, bool ReadOnly = false>
inline static IndexView<Value, NDims>
make_host_indexview(const Array& array) {
  typedef gridtools::storage_traits<BACKEND>::storage_info_t<0, NDims> storage_info_ty;
  typedef gridtools::storage_traits<BACKEND>::data_store_t<Value, storage_info_ty> data_store_t;

  data_store_t* ds = reinterpret_cast<data_store_t*>(const_cast<void*>(array.data()));

  return IndexView<Value, NDims>(gridtools::make_host_view(*ds));
}


// Old implementation
#else

template <typename Value, unsigned int NDims, bool ReadOnly = false>
inline static ArrayView<Value, NDims>
make_host_view(const Array& array) {
  return ArrayView<Value, NDims>((const Value*)(array.data()),array.shape());
}

template <typename Value, unsigned int NDims, bool ReadOnly = false>
static IndexView<Value, NDims>
make_host_indexview(const Array& array) {
  return IndexView<Value,NDims>( (Value*)(array.data()),array.shape().data() );
}

#endif

// --------------------------------------------------------------------------------------------

template <typename Value, unsigned int NDims, bool ReadOnly = false>
inline static IndexView<Value, NDims>
make_indexview(const Array& array) {
  impl::check_metadata<Value, NDims>(array);

  return make_host_indexview<Value,NDims>(array);
}

template <typename Value, unsigned int NDims, bool ReadOnly = false>
inline static ArrayView<Value, NDims>
make_view(const Array& array) {
    impl::check_metadata<Value, NDims>(array);

    return make_host_view<Value, NDims, ReadOnly>(array);
}

}
}

#endif
