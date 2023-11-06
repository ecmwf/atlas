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

#include "atlas/util/vector.h"

#include "atlas-io.h"

namespace atlas {

//---------------------------------------------------------------------------------------------------------------------

template <typename T>
void encode_data(const atlas::vector<T>& v, atlas::io::Data& out) {
    out.assign(v.data(), v.size() * sizeof(T));
}

//---------------------------------------------------------------------------------------------------------------------

template <typename T>
size_t encode_metadata(const atlas::vector<T>& v, atlas::io::Metadata& metadata) {
    using atlas::io::ArrayMetadata;
    using DataType = ArrayMetadata::DataType;
    return encode_metadata(ArrayMetadata{DataType::create<T>(), {v.size()}}, metadata);
}

//---------------------------------------------------------------------------------------------------------------------

template <typename T>
void decode(const atlas::io::Metadata& m, const atlas::io::Data& encoded, atlas::vector<T>& out) {
    atlas::io::ArrayMetadata array(m);
    const T* data = static_cast<const T*>(encoded.data());
    out.assign(data, data + array.size());
}

//---------------------------------------------------------------------------------------------------------------------

}  // end namespace atlas
