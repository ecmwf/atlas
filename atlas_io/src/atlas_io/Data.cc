/*
 * (C) Copyright 2020 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "Data.h"

#include <memory>

#include "eckit/utils/Compressor.h"

#include "atlas_io/Stream.h"
#include "atlas_io/Trace.h"
#include "atlas_io/detail/Assert.h"
#include "atlas_io/detail/Checksum.h"

namespace atlas {
namespace io {

//---------------------------------------------------------------------------------------------------------------------

Data::Data(void* p, size_t size): buffer_(p, size), size_(size) {}

std::uint64_t Data::write(Stream& out) const {
    ATLAS_IO_TRACE();
    if (size()) {
        ATLAS_IO_ASSERT(buffer_.size() >= size());
        return out.write(buffer_.data(), size());
    }
    return 0;
}

std::uint64_t Data::read(Stream& in, size_t size) {
    if (size > size_) {
        buffer_.resize(size);
        size_ = size;
    }
    return in.read(buffer_, size);
}


void Data::compress(const std::string& compression) {
    ATLAS_IO_TRACE("compress(" + compression + ")");
    if (size_) {
        auto compressor = std::unique_ptr<eckit::Compressor>(eckit::CompressorFactory::instance().build(compression));
        if (dynamic_cast<eckit::NoCompressor*>(compressor.get())) {
            return;
        }
        eckit::Buffer compressed(size_t(1.2 * size_));
        size_   = compressor->compress(buffer_, size_, compressed);
        buffer_ = std::move(compressed);
    }
}

void Data::decompress(const std::string& compression, size_t uncompressed_size) {
    ATLAS_IO_TRACE("decompress(" + compression + ")");

    auto compressor = std::unique_ptr<eckit::Compressor>(eckit::CompressorFactory::instance().build(compression));
    if (dynamic_cast<eckit::NoCompressor*>(compressor.get())) {
        return;
    }

    eckit::Buffer uncompressed(size_t(1.2 * uncompressed_size));
    compressor->uncompress(buffer_, size_, uncompressed, uncompressed_size);
    size_   = uncompressed_size;
    buffer_ = std::move(uncompressed);
}

void Data::clear() {
    buffer_ = eckit::Buffer{};
    size_   = 0;
}

std::string Data::checksum(const std::string& algorithm) const {
    return atlas::io::checksum(buffer_, size_, algorithm);
}

void Data::assign(const Data& other) {
    if (other.size() > buffer_.size()) {
        buffer_.resize(other.size());
    }
    size_ = other.size();
    buffer_.copy(other.buffer_, size_);
}

void Data::assign(const void* p, size_t s) {
    if (s > size()) {
        buffer_.resize(s);
    }
    size_ = s;
    buffer_.copy(p, size_);
}

Data& compress(Data& data, const std::string& compression) {
    data.compress(compression);
    return data;
}

Data& decompress(Data& data, const std::string& compression, size_t uncompressed_size) {
    data.decompress(compression, uncompressed_size);
    return data;
}

//---------------------------------------------------------------------------------------------------------------------

void encode(const Data& in, Data& out) {
    out.assign(in);
}

//---------------------------------------------------------------------------------------------------------------------

}  // namespace io
}  // namespace atlas
