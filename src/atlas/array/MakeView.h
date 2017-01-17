#pragma once

#include "atlas/array/array_fwd.h"
#include "atlas/array/ArrayView.h"
#include "atlas/array/IndexView.h"
#include "atlas/array/StorageView.h"
#include "atlas/internals/atlas_defines.h"

//------------------------------------------------------------------------------
#ifndef ATLAS_HAVE_GRIDTOOLS_STORAGE
// Old implementation

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

//------------------------------------------------------------------------------

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

// --------------------------------------------------------------------------------------------

} // namespace array
} // namespace atlas
#endif
