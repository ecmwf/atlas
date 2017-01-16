#pragma once

namespace atlas {
namespace array {

//#define GTNS gridtools::
//#define GTArray gridtools::Array
#define GTArray Array

class Array;

template<typename Value>
class StorageView;

template<typename Value, int RANK>
class ArrayView;

template<typename Value, int RANK>
class IndexView;


template <typename Value, unsigned int NDims, bool ReadOnly = false>
ArrayView<Value, NDims>
make_view(const Array& array);

template <typename Value, unsigned int NDims, bool ReadOnly = false>
ArrayView<Value, NDims>
make_host_view(const Array& array);

template <typename Value, unsigned int NDims, bool ReadOnly = false>
ArrayView<Value, NDims>
make_device_view(const Array& array);

template <typename Value, unsigned int NDims, bool ReadOnly = false>
IndexView<Value, NDims>
make_indexview(const Array& array);

template <typename Value, unsigned int NDims, bool ReadOnly = false>
IndexView<Value, NDims>
make_host_indexview(const Array& array);

template <typename Value>
StorageView<Value>
make_storageview(const Array& array);

template <typename Value>
StorageView<Value>
make_host_storageview(const Array& array);

template <typename Value>
StorageView<Value>
make_device_storageview(const Array& array);







//class ArrayBase;
//template <typename Value, unsigned int NDims, bool ReadOnly = false>
//ArrayView<Value, NDims>
//make_view(const ArrayBase& array);





}
}
