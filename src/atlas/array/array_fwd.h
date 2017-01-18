#pragma once

namespace atlas {
namespace array {

class DataType;

class ArraySpec;

class Array;

template <typename Value>
class ArrayT;

template <typename Value>
class StorageView;

template <typename Value, int RANK>
class ArrayView;

template <typename Value, int RANK>
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

}
}
