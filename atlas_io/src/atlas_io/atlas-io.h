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

#include <iomanip>
#include <sstream>
#include <string>

#include "atlas_io/detail/Link.h"
#include "atlas_io/detail/Reference.h"
#include "atlas_io/detail/StaticAssert.h"
#include "atlas_io/detail/sfinae.h"

#include "atlas_io/Exceptions.h"
#include "atlas_io/FileStream.h"
#include "atlas_io/Record.h"
#include "atlas_io/RecordItemReader.h"
#include "atlas_io/RecordPrinter.h"
#include "atlas_io/RecordReader.h"
#include "atlas_io/RecordWriter.h"
#include "atlas_io/Session.h"
#include "atlas_io/Stream.h"
#include "atlas_io/Trace.h"

#include "atlas_io/types/array.h"
#include "atlas_io/types/scalar.h"
#include "atlas_io/types/string.h"


namespace atlas {
namespace io {

//---------------------------------------------------------------------------------------------------------------------

inline Link link(const std::string& uri) {
    return Link{uri};
}

//---------------------------------------------------------------------------------------------------------------------

template <typename T, disable_if_interpretable_t<T, ArrayReference> = 0>
Reference<T> ref(const T& x, tag::enable_static_assert = tag::enable_static_assert()) {
    static_assert(is_encodable<T>(),
                  "in atlas::io::ref(const Value&)"
                  "\n"
                  "\n     Static assertion failed"
                  "\n     -----------------------"
                  "\n"
                  "\n     Cannot encode values of referenced type."
                  "\n"
                  "\n     Implement the functions"
                  "\n"
                  "\n         void encode_data(const Value& in, atlas::io::Data& out);"
                  "\n         size_t encode_metadata(const Value& value, atlas::io::Metadata& metadata);"
                  "\n"
                  "\n     or alternatively a conversion function to atlas::io::types::ArrayView"
                  "\n"
                  "\n         void interprete(const Value& in, atlas::io::types::ArrayView& out)"
                  "\n"
                  "\n     Rules of argument-dependent-lookup apply."
                  "\n     --> Functions need to be declared in namespace of any of the arguments."
                  "\n"
                  "\n     Note, turn this into a runtime exception by calling this function instead:"
                  "\n"
                  "\n        atlas::io::ref(const T&, atlas::io::no_static_assert() )"
                  "\n");
    return Reference<T>(x);
}


template <typename T, disable_if_interpretable_t<T, ArrayReference> = 0>
Reference<T> ref(const T& x, tag::disable_static_assert) {
    if (not is_encodable<T>()) {
        throw NotEncodable(x);
    }
    return Reference<T>(x);
}


template <typename T, enable_if_interpretable_t<T, ArrayReference> = 0>
ArrayReference ref(const T& x, tag::enable_static_assert = tag::enable_static_assert()) {
    ArrayReference w;
    interprete(x, w);
    return w;
}

//---------------------------------------------------------------------------------------------------------------------

template <typename T>
RecordItem copy(T&& value, tag::disable_static_assert) {
    return RecordItem(std::forward<T>(value), tag::disable_static_assert());
}

template <typename T>
RecordItem copy(T&& value) {
    return RecordItem(std::forward<T>(value));
}

//---------------------------------------------------------------------------------------------------------------------

template <typename T>
void encode(const T& in, atlas::io::Metadata& metadata, atlas::io::Data& data,
            tag::enable_static_assert = tag::enable_static_assert()) {
    auto referenced = ref(in, tag::enable_static_assert());
    sfinae::encode_metadata(referenced, metadata);
    sfinae::encode_data(referenced, data);
}

template <typename T>
void encode(const T& in, atlas::io::Metadata& metadata, atlas::io::Data& data, tag::disable_static_assert) {
    auto referenced = ref(in, tag::disable_static_assert());
    sfinae::encode_metadata(referenced, metadata);
    sfinae::encode_data(referenced, data);
}

//---------------------------------------------------------------------------------------------------------------------

}  // namespace io
}  // namespace atlas
