/*
 * (C) Copyright 2020 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#pragma once

#include "atlas_io/detail/TypeTraits.h"

namespace atlas {
namespace io {

// -------------------------------------------------------------------------------------------------------

namespace {

// These anonymous namespace functions are to avoid recursive behaviour in
// following sfinae namespace

template <typename T, typename A, enable_if_interpretable_t<T, A> = 0>
inline void do_interprete(const T& in, A& interpreted) {
    interprete(in, interpreted);
}

template <typename T, enable_if_can_encode_metadata_t<T> = 0>
inline size_t do_encode_metadata(const T& in, Metadata& out) {
    size_t size = encode_metadata(in, out);
    return size;
}

template <typename T, enable_if_can_encode_data_t<T> = 0>
inline void do_encode_data(const T& in, Data& out) {
    encode_data(in, out);
}

}  // namespace

// -------------------------------------------------------------------------------------------------------

namespace sfinae {

// -------------------------------------------------------------------------------------------------------

template <typename T, typename A, enable_if_interpretable_t<T, A> = 0>
bool interprete(const T& in, A& interpreted) {
    do_interprete(in, interpreted);
    return true;
}

template <typename T, typename A, disable_if_interpretable_t<T, A> = 0>
bool interprete(const T& /*in*/, A& /*interpreted*/) {
    return false;
}

// -------------------------------------------------------------------------------------------------------

template <typename T, enable_if_can_encode_metadata_t<T> = 0>
bool encode_data(const T& in, Data& out) {
    do_encode_data(in, out);
    return true;
}

template <typename T, disable_if_can_encode_metadata_t<T> = 0>
bool encode_data(const T&, Data&) {
    return false;
}

// -------------------------------------------------------------------------------------------------------

template <typename T, enable_if_can_encode_data_t<T> = 0>
bool encode_metadata(const T& in, Metadata& out) {
    do_encode_metadata(in, out);
    return true;
}

template <typename T, disable_if_can_encode_data_t<T> = 0>
bool encode_metadata(const T&, Metadata&) {
    return false;
}

template <typename T, enable_if_can_encode_data_t<T> = 0>
bool encode_metadata(const T& in, Metadata& out, size_t& data_size) {
    data_size = do_encode_metadata(in, out);
    return true;
}

template <typename T, disable_if_can_encode_data_t<T> = 0>
bool encode_metadata(const T&, Metadata&, size_t& data_size) {
    data_size = 0;
    return false;
}

// -------------------------------------------------------------------------------------------------------

}  // namespace sfinae

}  // namespace io
}  // namespace atlas
