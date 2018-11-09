/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

/// @author Willem Deconinck

#pragma once

#include <cstddef>
#include "atlas/array/ArrayViewDefs.h"


namespace atlas {
namespace array {

class DataType;

class ArraySpec;

class ArrayShape;

class ArrayStrides;

class Array;

template <typename Value>
class ArrayT;

template <typename Value, int RANK, Intent AccessMode>
class ArrayView;

template <typename Value, int RANK, Intent AccessMode>
class LocalView;

template <typename Value, int RANK>
class IndexView;

template <typename Value, unsigned int NDims, Intent AccessMode = Intent::ReadWrite>
ArrayView<Value, NDims, AccessMode> make_view( const Array& array );

template <typename Value, unsigned int NDims, Intent AccessMode = Intent::ReadWrite>
LocalView<Value, NDims, AccessMode> make_view( const Value data[], const ArrayShape& );

template <typename Value, unsigned int NDims, Intent AccessMode = Intent::ReadWrite>
LocalView<Value, NDims, AccessMode> make_view( const Value data[], size_t );

template <typename Value, unsigned int NDims, Intent AccessMode = Intent::ReadWrite>
ArrayView<Value, NDims, AccessMode> make_host_view( const Array& array );

template <typename Value, unsigned int NDims, Intent AccessMode = Intent::ReadWrite>
ArrayView<Value, NDims, AccessMode> make_device_view( const Array& array );

template <typename Value, unsigned int NDims, Intent AccessMode = Intent::ReadWrite>
IndexView<Value, NDims> make_indexview( const Array& array );

template <typename Value, unsigned int NDims, Intent AccessMode = Intent::ReadWrite>
IndexView<Value, NDims> make_host_indexview( const Array& array );

class Table;

template <bool ReadOnly>
class TableView;

template <bool ReadOnly>
class TableRow;

template <bool ReadOnly = true>
TableView<ReadOnly> make_table_view( const Table& table );

}  // namespace array
}  // namespace atlas
