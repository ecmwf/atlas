/*
 * (C) Copyright 2024- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "trace.h"

#include <cstdlib>
#include <sstream>
#include <iomanip>

namespace pluto::trace {

static bool PLUTO_TRACE() {
    char* val;
    val = std::getenv("PLUTO_TRACE");
    if (val != nullptr) {
        return std::atoi(val);
    }
    return false;
}

static bool PLUTO_TRACE_FORMAT_BYTES() {
    char* val;
    val = std::getenv("PLUTO_TRACE_FORMAT_BYTES");
    if (val != nullptr) {
        return std::atoi(val);
    }
    return true;
}

Options& options() {
    static Options opts{PLUTO_TRACE()};
    return opts;
}

OutputStream out;

std::string format_bytes(std::size_t bytes) {
    constexpr double KB = 1024.;
    constexpr double MB = 1024. * KB;
    constexpr double GB = 1024. * MB;

    std::stringstream ss;
    if (PLUTO_TRACE_FORMAT_BYTES()) {
        ss << std::setprecision(2) << std::fixed;
        double b = static_cast<double>(bytes);
        if (b >= GB) {
            ss << b / GB << "G";
        }
        else if (b >= MB) {
            ss << b / MB << "M";
        }
        else if (b >= KB) {
            ss << b / KB << "K";
        }
        else {
            ss << b << "B";
        }
    }
    else {
        ss << bytes;
    }
    return ss.str();
}

}  // namespace pluto::trace
