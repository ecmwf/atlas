/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

/// @author Willem Deconinck

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

template <typename Value, int RANK, bool ReadOnly = false>
class ArrayView;

template <typename Value, int RANK>
class IndexView;

template <typename Value, unsigned int NDims, bool ReadOnly = false>
ArrayView<Value, NDims, ReadOnly>
make_view(const Array& array);

template <typename Value, unsigned int NDims, bool ReadOnly = false>
ArrayView<Value, NDims, ReadOnly>
make_host_view(const Array& array);

template <typename Value, unsigned int NDims, bool ReadOnly = false>
ArrayView<Value, NDims, ReadOnly>
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

class Table;

template <bool ReadOnly>
class TableView;

template <bool ReadOnly>
class TableRow;

template<bool ReadOnly=true>
TableView<ReadOnly>
make_table_view(const Table& table);

}
}
