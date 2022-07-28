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

namespace eckit {
class JSON;
}
namespace atlas {
namespace io {

/// Store UTC time up to nanosecond precision
struct Time {
    std::uint64_t tv_sec{0};   ///<    seconds since Epoch (1970-01-01T00:00:00Z)
    std::uint64_t tv_nsec{0};  ///<    additional nanoseconds

    /// Create current time using system clock
    static Time now();

    /// print UTC time in ISO 8601 format: "1970-01-01T00:00:00.123456789Z"
    void print(std::ostream&) const;

    friend std::ostream& operator<<(std::ostream&, const Time&);
    friend eckit::JSON& operator<<(eckit::JSON&, const Time&);

    /// @return string of UTC time in ISO 8601 format: "1970-01-01T00:00:00.123456789Z"
    std::string str() const;
};


}  // namespace io
}  // namespace atlas
