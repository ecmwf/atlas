/*
 * (C) Copyright 2023- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#pragma once

#include "eckit/codec/codec.h"

namespace atlas::io {

    // Encoding/Decoding
    using Metadata = ::eckit::codec::Metadata;
    using DataType = ::eckit::codec::DataType;
    using Data     = ::eckit::codec::Data;
    using Encoder  = ::eckit::codec::Encoder;
    using Decoder  = ::eckit::codec::Decoder;

    // Record
    using Record           = ::eckit::codec::Record;
    using RecordWriter     = ::eckit::codec::RecordWriter;
    using RecordPrinter    = ::eckit::codec::RecordPrinter;
    using RecordItemReader = ::eckit::codec::RecordItemReader;
    using RecordReader     = ::eckit::codec::RecordReader;

    // I/O
    using Session          = ::eckit::codec::Session;
    using Mode             = ::eckit::codec::Mode;
    using Stream           = ::eckit::codec::Stream;
    using FileStream       = ::eckit::codec::FileStream;
    using InputFileStream  = ::eckit::codec::InputFileStream;
    using OutputFileStream = ::eckit::codec::OutputFileStream;

    // Array
    using ArrayReference = ::eckit::codec::ArrayReference;
    using ArrayMetadata  = ::eckit::codec::ArrayMetadata;
    using ArrayShape     = ::eckit::codec::ArrayShape;

    // Exceptions
    using Exception      = ::eckit::codec::Exception;
    using NotEncodable   = ::eckit::codec::NotEncodable;
    using NotDecodable   = ::eckit::codec::NotDecodable;
    using InvalidRecord  = ::eckit::codec::InvalidRecord;
    using DataCorruption = ::eckit::codec::DataCorruption;
    using WriteError     = ::eckit::codec::WriteError;

    // Type traits
    template <class T, class A>
    static constexpr bool is_interpretable() {
        return ::eckit::codec::is_interpretable<T,A>();
    }
    template <typename T>
    static constexpr bool is_encodable() {
        return ::eckit::codec::is_encodable<T>();
    }
    template <typename T>
    static constexpr bool is_decodable() {
        return ::eckit::codec::is_decodable<T>();
    }
    template <typename T>
    static constexpr bool can_encode_metadata() {
        return ::eckit::codec::can_encode_metadata<T>();
    }
    template <typename T>
    static constexpr bool can_encode_data() {
        return ::eckit::codec::can_encode_metadata<T>();
    }

    namespace tag {
        using disable_static_assert = ::eckit::codec::tag::disable_static_assert;
        using enable_static_assert  = ::eckit::codec::tag::enable_static_assert;
    }

    // Functions
    using ::eckit::codec::ref;
    using ::eckit::codec::copy;
    using ::eckit::codec::encode;
    using ::eckit::codec::decode;
    using ::eckit::codec::interprete;
    using ::eckit::codec::link;
    using ::eckit::codec::make_datatype;
    // template <typename DATATYPE>
    // using make_datatype = eckit::codec::make_datatype<DATATYPE>;

    namespace defaults {
        using ::eckit::codec::defaults::compression_algorithm;
        using ::eckit::codec::defaults::checksum_algorithm;
        using ::eckit::codec::defaults::checksum_read;
        using ::eckit::codec::defaults::checksum_write;
    }
}

#define ATLAS_IO_ASSERT(X) ASSERT(X)
#define ATLAS_IO_ASSERT_MSG(X, M) ASSERT_MSG(X, M)
