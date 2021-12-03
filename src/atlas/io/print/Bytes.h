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
#include <iosfwd>
#include <string>

namespace atlas {
namespace io {

struct Bytes {
public:
    Bytes(size_t bytes): bytes_(bytes) {}

    operator size_t() const { return bytes_; }

    std::string str(int decimals = 2, int width = 7) const;

    void print(std::ostream& out, int decimals = 2, int width = 7) const;

    friend std::ostream& operator<<(std::ostream&, const Bytes&);

private:
    size_t bytes_;
};

}  // namespace io
}  // namespace atlas
