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

#include <string>

#include "eckit/config/Resource.h"

namespace atlas {
namespace io {
namespace defaults {

[[maybe_unused]] static std::string checksum_algorithm() {
    static std::string checksum =
        eckit::Resource<std::string>("atlas.io.checksum.algorithm;$ATLAS_IO_CHECKSUM", "xxh64");
    return checksum;
}

[[maybe_unused]] static bool checksum_read() {
    static bool checksum = eckit::Resource<bool>("atlas.io.checksum.read;$ATLAS_IO_CHECKSUM_READ", true);
    return checksum;
}

[[maybe_unused]] static bool checksum_write() {
    static bool checksum = eckit::Resource<bool>("atlas.io.checksum.write;$ATLAS_IO_CHECKSUM_WRITE", true);
    return checksum;
}

[[maybe_unused]] static const std::string& compression_algorithm() {
    static std::string compression = eckit::Resource<std::string>("atlas.io.compression;$ATLAS_IO_COMPRESSION", "none");
    return compression;
}


}  // namespace defaults
}  // namespace io
}  // namespace atlas
