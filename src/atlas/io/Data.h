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

#include "eckit/io/Buffer.h"

namespace atlas {
namespace io {

class Stream;

//---------------------------------------------------------------------------------------------------------------------

class Data {
public:
    Data() = default;
    Data(void*, size_t);
    Data(Data&&)  = default;
    Data& operator=(Data&&) = default;

    operator const void*() const { return data(); }
    const void* data() const { return buffer_.data(); }
    size_t size() const { return size_; }

    void assign(const Data& other);
    void assign(const void*, size_t);
    void clear();

    std::uint64_t write(Stream& out) const;
    std::uint64_t read(Stream& in, size_t size);
    void compress(const std::string& compression);
    void decompress(const std::string& compression, size_t uncompressed_size);
    std::string checksum(const std::string& algorithm = "") const;

private:
    eckit::Buffer buffer_;
    size_t size_{0};
};

//---------------------------------------------------------------------------------------------------------------------

void encode_data(const Data&, Data& out);

//---------------------------------------------------------------------------------------------------------------------

}  // namespace io
}  // namespace atlas
