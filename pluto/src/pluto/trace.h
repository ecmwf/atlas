/*
 * (C) Copyright 2024- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */
#pragma once

#include <cstdlib>
#include <iostream>

namespace pluto {
struct TraceOptions {
    static TraceOptions& instance() {
        static TraceOptions opts;
        return opts;
    }
    bool enabled{false};
    std::ostream* out{&std::cout};

private:
    TraceOptions() {
        char* val;
        val = std::getenv("PLUTO_TRACE");
        if (val != NULL) {
            enabled = std::atoi(val);
        }
    }
};

inline std::ostream& trace() {
    return *TraceOptions::instance().out;
}

inline bool trace_enabled() {
    return TraceOptions::instance().enabled;
}

inline void set_trace(bool value) {
    TraceOptions::instance().enabled = value;
}

}  // namespace pluto
