/*
 * (C) Copyright 2020 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

// Cray C++ compiler should not try to optimize this code
#ifdef _CRAYC
#pragma _CRI noopt
#endif

// GNU C++ compiler (version 11) should not try to optimize this code
#ifdef __GNUC__
#pragma GCC optimize("O0")
#endif

#include "scalar.h"

#include <cstdint>
#include <string>

#include "eckit/utils/ByteSwap.h"

#include "atlas_io/atlas_compat.h"
#include "atlas_io/detail/Assert.h"
#include "atlas_io/detail/DataType.h"

#include "atlas_io/Data.h"
#include "atlas_io/Metadata.h"
#include "atlas_io/detail/Base64.h"

namespace atlas {
namespace io {


//---------------------------------------------------------------------------------------------------------------------

template <typename T>
void decode_scalar(const atlas::io::Metadata& metadata, T& value) {
    ATLAS_IO_ASSERT(metadata.getString("type") == "scalar");
    ATLAS_IO_ASSERT(metadata.getString("datatype") == DataType::str<T>());
    metadata.get("value", value);
}

template <typename T>
void decode_scalar_b64(const atlas::io::Metadata& metadata, T& value) {
    ATLAS_IO_ASSERT(metadata.getString("type") == "scalar");
    ATLAS_IO_ASSERT(metadata.getString("datatype") == DataType::str<T>());
    std::string base64 = metadata.getString("base64");
    T value_ns         = Base64::decode<T>(base64);
    if (Endian::native == Endian::little) {
        T tmp = value_ns;
        eckit::byteswap(tmp);
        value_ns = tmp;
    }
    value = value_ns;
}

//---------------------------------------------------------------------------------------------------------------------

template <typename T>
void encode_scalar_metadata(const T& value, atlas::io::Metadata& out) {
    out.set("type", "scalar");
    out.set("datatype", DataType::str<T>());
    out.set("value", value);
}

inline void encode_scalar_metadata(const unsigned long& value, atlas::io::Metadata& out) {
    out.set("type", "scalar");
    out.set("datatype", DataType::str<size_t>());
    out.set("value", size_t(value));
}

inline void encode_scalar_metadata(const unsigned long long& value, atlas::io::Metadata& out) {
    out.set("type", "scalar");
    out.set("datatype", DataType::str<size_t>());
    out.set("value", size_t(value));
}

template <typename T>
size_t encode_scalar_metadata_b64(const T& value, atlas::io::Metadata& out) {
    encode_scalar_metadata(value, out);
    T value_ns = value;
    if (Endian::native == Endian::little) {
        eckit::byteswap(value_ns);
    }
    out.set("base64", Base64::encode(value_ns));
    return 0;
}

//---------------------------------------------------------------------------------------------------------------------

size_t encode_metadata(const int& value, atlas::io::Metadata& out) {
    return encode_scalar_metadata_b64(value, out);
}

size_t encode_metadata(const long& value, atlas::io::Metadata& out) {
    return encode_scalar_metadata_b64(value, out);
}

size_t encode_metadata(const long long& value, atlas::io::Metadata& out) {
    return encode_scalar_metadata_b64(value, out);
}

size_t encode_metadata(const unsigned long& value, atlas::io::Metadata& out) {
    return encode_scalar_metadata_b64(value, out);
}

size_t encode_metadata(const unsigned long long& value, atlas::io::Metadata& out) {
    return encode_scalar_metadata_b64(value, out);
}

size_t encode_metadata(const float& value, atlas::io::Metadata& out) {
    return encode_scalar_metadata_b64(value, out);
}

size_t encode_metadata(const double& value, atlas::io::Metadata& out) {
    return encode_scalar_metadata_b64(value, out);
}

//---------------------------------------------------------------------------------------------------------------------

void encode_data(const int&, atlas::io::Data&) {}
void encode_data(const long&, atlas::io::Data&) {}
void encode_data(const long long&, atlas::io::Data&) {}
void encode_data(const unsigned long&, atlas::io::Data&) {}
void encode_data(const unsigned long long&, atlas::io::Data&) {}
void encode_data(const float&, atlas::io::Data&) {}
void encode_data(const double&, atlas::io::Data&) {}

//---------------------------------------------------------------------------------------------------------------------

void decode(const atlas::io::Metadata& metadata, const atlas::io::Data&, int& value) {
    decode_scalar_b64(metadata, value);
}
void decode(const atlas::io::Metadata& metadata, const atlas::io::Data&, long& value) {
    decode_scalar_b64(metadata, value);
}
void decode(const atlas::io::Metadata& metadata, const atlas::io::Data&, long long& value) {
    decode_scalar_b64(metadata, value);
}
void decode(const atlas::io::Metadata& metadata, const atlas::io::Data&, unsigned long& value) {
    decode_scalar_b64(metadata, value);
}
void decode(const atlas::io::Metadata& metadata, const atlas::io::Data&, unsigned long long& value) {
    decode_scalar_b64(metadata, value);
}
void decode(const atlas::io::Metadata& metadata, const atlas::io::Data&, float& value) {
    decode_scalar_b64(metadata, value);
}
void decode(const atlas::io::Metadata& metadata, const atlas::io::Data&, double& value) {
    decode_scalar_b64(metadata, value);
}

//---------------------------------------------------------------------------------------------------------------------

}  // namespace io
}  // namespace atlas
