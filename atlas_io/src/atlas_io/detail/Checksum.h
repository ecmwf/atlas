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

namespace atlas {
namespace io {

class Checksum {
public:
    Checksum() = default;
    Checksum(const std::string& checksum);
    bool available() const;
    std::string str() const;
    std::string str(size_t size) const;
    std::string algorithm() const { return algorithm_; }

private:
    std::string algorithm_;
    std::string checksum_;
};

std::string checksum(const void* buffer, size_t size, const std::string& algorithm = "");


}  // namespace io
}  // namespace atlas
