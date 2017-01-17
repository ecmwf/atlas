#pragma once

#include "atlas/array/array_fwd.h"
#include "atlas/array/ArrayView.h"
#include "atlas/array/IndexView.h"
#include "atlas/array/StorageView.h"
#include "atlas/internals/atlas_defines.h"

//------------------------------------------------------------------------------
#ifndef ATLAS_HAVE_GRIDTOOLS_STORAGE
#include "atlas/array/Array.h"

// Native implementation

namespace atlas {
namespace array {

namespace {
    template<typename Value, unsigned Rank>
    inline static void check_metadata(const Array& array)
    {
        if(array.rank() != Rank ) {
            std::stringstream err;
            err << "Number of dimensions do not match: template argument " << Rank << " expected to be " << array.rank();
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

template <typename Value, unsigned int Rank, bool ReadOnly>
inline ArrayView<Value, Rank>
make_host_view(const Array& array) {
    return ArrayView<Value, Rank>((const Value*)(array.storage()),array.shape());
}


template <typename Value, unsigned int Rank, bool ReadOnly>
inline ArrayView<Value, Rank>
make_device_view(const Array& array) {
    return make_host_view<Value,Rank,ReadOnly>(array);
}


template <typename Value, unsigned int Rank, bool ReadOnly>
inline IndexView<Value, Rank>
make_host_indexview(const Array& array) {
    return IndexView<Value,Rank>( (Value*)(array.storage()),array.shape().data() );
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


template <typename Value, unsigned int Rank, bool ReadOnly>
inline IndexView<Value, Rank>
make_indexview(const Array& array) {
    check_metadata<Value, Rank>(array);
    return make_host_indexview<Value,Rank>(array);
}


template <typename Value, unsigned int Rank, bool ReadOnly>
inline ArrayView<Value, Rank>
make_view(const Array& array) {
    check_metadata<Value, Rank>(array);
    return make_host_view<Value, Rank, ReadOnly>(array);
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
