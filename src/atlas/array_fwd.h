/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

/// @file array_fwd.h
/// @brief Forward declarations of public API header atlas/array.h
/// @author Willem Deconinck

/// @file array_fwd.h
/// @author Willem Deconinck

#pragma once

#include <cstddef>
#include <type_traits>
#include "atlas/array/ArrayViewDefs.h"

#define ENABLE_IF_NOT_CONST typename std::enable_if<!std::is_const<Value>::value, Value>::type* = nullptr
#define ENABLE_IF_CONST typename std::enable_if<std::is_const<Value>::value, Value>::type* = nullptr

namespace atlas {
namespace array {

class DataType;

class ArraySpec;

class ArrayShape;

class ArrayStrides;

class Array;

template <typename Value>
class ArrayT;

template <typename Value, int Rank>
class ArrayView;

template <typename Value, int Rank>
class LocalView;

template <typename Value, int Rank>
class IndexView;

template <typename Value, int Rank>
ArrayView<Value, Rank> make_view(Array& array);

template <typename Value, int Rank>
ArrayView<const Value, Rank> make_view(const Array& array);


template <class Value, int Rank, ENABLE_IF_NOT_CONST>
LocalView<Value, Rank> make_view(Value data[], const ArrayShape& shape);

template <class Value, int Rank, ENABLE_IF_NOT_CONST>
LocalView<const Value, Rank> make_view(const Value data[], const ArrayShape& shape);

template <class Value, int Rank, ENABLE_IF_CONST>
LocalView<Value, Rank> make_view(Value data[], const ArrayShape& shape);

template <class Value, int Rank, ENABLE_IF_CONST>
LocalView<Value, Rank> make_view(typename std::remove_const<Value>::type data[], const ArrayShape& shape);

//------------------------------------------------------------------------------------------------------

template <typename Value, unsigned int Rank, ENABLE_IF_NOT_CONST>
LocalView<Value, Rank> make_view(Value data[], size_t size);

template <typename Value, unsigned int Rank, ENABLE_IF_NOT_CONST>
LocalView<const Value, Rank> make_view(const Value data[], size_t size);

template <typename Value, unsigned int Rank, ENABLE_IF_CONST>
LocalView<Value, Rank> make_view(Value data[], size_t size);

template <typename Value, unsigned int Rank, ENABLE_IF_CONST>
LocalView<Value, Rank> make_view(typename std::remove_const<Value>::type data[], size_t size);

template <typename Value, int Rank>
ArrayView<Value, Rank> make_host_view(Array& array);

template <typename Value, int Rank>
ArrayView<const Value, Rank> make_host_view(const Array& array);

template <typename Value, int Rank>
ArrayView<Value, Rank> make_device_view(Array& array);

template <typename Value, int Rank>
ArrayView<const Value, Rank> make_device_view(const Array& array);


template <typename Value, int Rank>
IndexView<Value, Rank> make_indexview(Array& array);

template <typename Value, int Rank>
IndexView<const Value, Rank> make_indexview(const Array& array);


template <typename Value, int Rank>
IndexView<Value, Rank> make_host_indexview(Array& array);

template <typename Value, int Rank>
IndexView<const Value, Rank> make_host_indexview(const Array& array);

// class Table;

// template <bool ReadOnly>
// class TableView;

// template <bool ReadOnly>
// class TableRow;

// template <bool ReadOnly = true>
// TableView<ReadOnly> make_table_view( const Table& table );

#undef ENABLE_IF_NOT_CONST
#undef ENABLE_IF_CONST


}  // namespace array
}  // namespace atlas
