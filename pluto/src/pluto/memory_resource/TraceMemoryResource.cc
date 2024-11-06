/*
 * (C) Copyright 2024- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "TraceMemoryResource.h"

#include <limits>
#include <sstream>
#include <iomanip>

namespace pluto {

static std::string bytes_to_string(std::size_t bytes) {
    constexpr double KB = 1024.;
    constexpr double MB = 1024.*KB;
    constexpr double GB = 1024.*MB;
    double b = static_cast<double>(bytes);

    std::stringstream ss;
    ss << std::setprecision(2) << std::fixed;
    if (b >= GB) {
        ss << b/GB << "G";
    }
    else if (b >= MB) {
        ss << b/MB << "M";
    }
    else if (b >= KB) {
        ss << b/KB << "K";
    }
    else {
        ss << b << "B";
    }
    return ss.str();
}


int TraceMemoryResource::nest = 0;

void* TraceMemoryResource::do_allocate(std::size_t bytes, std::size_t alignment) {
    nest++;
    if (TraceOptions::instance().enabled) {
        *TraceOptions::instance().out << std::string(4*nest,' ') << "[" << name_ << " (alloc)] " << bytes << " bytes (" << bytes_to_string(bytes) << ") with alignment " << alignment << std::endl;
    }
    auto ptr = mr_->allocate(bytes, alignment);
    nest--;
    return  ptr;
}

void TraceMemoryResource::do_deallocate(void* p, std::size_t bytes, std::size_t alignment) {
    nest++;
    // print_stacktrace(out_, 4*nest, CURRENT_STACKTRACE());
    if (TraceOptions::instance().enabled) {
        if (bytes == std::numeric_limits<std::size_t>::max()) {
            *TraceOptions::instance().out << std::string(4*nest,' ') << "[" << name_ << " (dealloc)] " << p << std::endl;
        }
        else {
            *TraceOptions::instance().out << std::string(4*nest,' ') << "[" << name_ << " (dealloc)] " << p << " : " << bytes << " bytes (" << bytes_to_string(bytes) << ") with alignment " << alignment << std::endl;
        }
    }
    mr_->deallocate(p, bytes, alignment);
    nest--;
}

}
