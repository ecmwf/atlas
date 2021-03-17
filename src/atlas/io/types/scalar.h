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

#include <cstdint>
#include <string>

#include "eckit/utils/ByteSwap.h"

#include "atlas/array/DataType.h"
#include "atlas/runtime/Exception.h"

#include "atlas/io/Data.h"
#include "atlas/io/Metadata.h"
#include "atlas/io/detail/Base64.h"

namespace atlas {
namespace io {

//---------------------------------------------------------------------------------------------------------------------

template <typename T>
void decode_scalar( const atlas::io::Metadata& metadata, T& value ) {
    ATLAS_ASSERT( metadata.getString( "type" ) == "scalar" );
    ATLAS_ASSERT( metadata.getString( "datatype" ) == array::DataType::str<T>() );
    metadata.get( "value", value );
}

template <typename T>
void decode_scalar_b64( const atlas::io::Metadata& metadata, T& value ) {
    ATLAS_ASSERT( metadata.getString( "type" ) == "scalar" );
    ATLAS_ASSERT( metadata.getString( "datatype" ) == array::DataType::str<T>() );
    std::string base64 = metadata.getString( "base64" );
    T value_ns         = Base64::decode<T>( base64 );
    if ( Endian::native == Endian::little ) {
        eckit::byteswap( value_ns );
    }
    value = value_ns;
}

//---------------------------------------------------------------------------------------------------------------------

template <typename T>
size_t encode_scalar_metadata( const T& value, atlas::io::Metadata& out ) {
    out.set( "type", "scalar" );
    out.set( "datatype", array::DataType::str<T>() );
    out.set( "value", value );
    return 0;
}

template <typename T>
size_t encode_scalar_metadata_b64( const T& value, atlas::io::Metadata& out ) {
    encode_scalar_metadata( value, out );
    T value_ns = value;
    if ( Endian::native == Endian::little ) {
        eckit::byteswap( value_ns );
    }
    out.set( "base64", Base64::encode( value_ns ) );
    return 0;
}

//---------------------------------------------------------------------------------------------------------------------

size_t encode_metadata( const std::int32_t& value, atlas::io::Metadata& out ) {
    return encode_scalar_metadata_b64( value, out );
}

size_t encode_metadata( const std::int64_t& value, atlas::io::Metadata& out ) {
    return encode_scalar_metadata_b64( value, out );
}

size_t encode_metadata( const size_t& value, atlas::io::Metadata& out ) {
    return encode_scalar_metadata_b64( value, out );
}

size_t encode_metadata( const float& value, atlas::io::Metadata& out ) {
    return encode_scalar_metadata_b64( value, out );
}

size_t encode_metadata( const double& value, atlas::io::Metadata& out ) {
    return encode_scalar_metadata_b64( value, out );
}

//---------------------------------------------------------------------------------------------------------------------

void encode_data( const std::int32_t&, atlas::io::Data& ) {}
void encode_data( const std::int64_t&, atlas::io::Data& ) {}
void encode_data( const float&, atlas::io::Data& ) {}
void encode_data( const double&, atlas::io::Data& ) {}
void encode_data( const size_t&, atlas::io::Data& ) {}

//---------------------------------------------------------------------------------------------------------------------

void decode( const atlas::io::Metadata& metadata, const atlas::io::Data&, std::int32_t& value ) {
    decode_scalar_b64( metadata, value );
}
void decode( const atlas::io::Metadata& metadata, const atlas::io::Data&, std::int64_t& value ) {
    decode_scalar_b64( metadata, value );
}
void decode( const atlas::io::Metadata& metadata, const atlas::io::Data&, size_t& value ) {
    decode_scalar_b64( metadata, value );
}
void decode( const atlas::io::Metadata& metadata, const atlas::io::Data&, float& value ) {
    decode_scalar_b64( metadata, value );
}
void decode( const atlas::io::Metadata& metadata, const atlas::io::Data&, double& value ) {
    decode_scalar_b64( metadata, value );
}

//---------------------------------------------------------------------------------------------------------------------

}  // namespace io
}  // namespace atlas
