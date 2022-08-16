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


#include "atlas_io/detail/Checksum.h"
#include "atlas_io/detail/Endian.h"

namespace atlas {
namespace io {

class DataInfo {
public:
    int section() const { return section_; }
    const std::string& compression() const { return compression_; }
    Endian endian() const { return endian_; }

    void endian(Endian e) { endian_ = e; }

    void compression(const std::string& c) { compression_ = c; }
    void size(size_t s) { uncompressed_size_ = s; }
    size_t size() const { return uncompressed_size_; }
    void compressed_size(size_t s) { compressed_size_ = s; }
    size_t compressed_size() const { return compressed_size_; }
    void compressed(bool f) {
        if (f == false) {
            compression("none");
        }
    }

    bool compressed() const { return compression_ != "none"; }

    operator bool() const { return section_ > 0; }

    const Checksum& checksum() const { return checksum_; }
    void checksum(const std::string& c) { checksum_ = Checksum(c); }

    void section(int s) { section_ = s; }

private:
    int section_{0};
    std::string compression_{"none"};
    Checksum checksum_;
    Endian endian_{Endian::native};
    size_t uncompressed_size_{0};
    size_t compressed_size_{0};
};

}  // namespace io
}  // namespace atlas
