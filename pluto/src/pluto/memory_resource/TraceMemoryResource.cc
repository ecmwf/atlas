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

#include "pluto/offload/Stream.h"

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
        *TraceOptions::instance().out << std::string(4*nest,' ') << "[" << name_ << " (alloc)] { bytes:" << bytes_to_string(bytes) << ", alignment:" << alignment << " }" << std::endl;
    }
    auto ptr = mr_->allocate(bytes, alignment);
    nest--;
    return  ptr;
}

void TraceMemoryResource::do_deallocate(void* p, std::size_t bytes, std::size_t alignment) {
    nest++;
    if (TraceOptions::instance().enabled) {
        if (bytes == std::numeric_limits<std::size_t>::max()) {
            *TraceOptions::instance().out << std::string(4*nest,' ') << "[" << name_ << " (dealloc)] { pointer:" << p << " }" << std::endl;
        }
        else {
            *TraceOptions::instance().out << std::string(4*nest,' ') << "[" << name_ << " (dealloc)] { pointer:" << p << ", bytes:" << bytes_to_string(bytes) << ", alignment:" << alignment << " }" << std::endl;
        }
    }
    mr_->deallocate(p, bytes, alignment);
    nest--;
}

void* TraceMemoryResource::do_allocate_async(std::size_t bytes, std::size_t alignment, const Stream& stream) {
    nest++;
    auto* async_mr = dynamic_cast<async_memory_resource*>(mr_);

    if (TraceOptions::instance().enabled) {
        *TraceOptions::instance().out << std::string(4*nest,' ') << "[" << name_ << " (alloc_async)] { bytes:" << bytes_to_string(bytes) << ", alignment:" << alignment << ", stream:" << stream.value() << " }" << std::endl;
    }
    void* ptr;
    if (async_mr) {
        ptr = async_mr->allocate_async(bytes, alignment, stream);
    }
    else {
        ptr = mr_->allocate(bytes, alignment);
    }
    nest--;
    return  ptr;
}

void TraceMemoryResource::do_deallocate_async(void* p, std::size_t bytes, std::size_t alignment, const Stream& stream) {
    nest++;
    auto* async_mr = dynamic_cast<async_memory_resource*>(mr_);
    if (TraceOptions::instance().enabled) {
        if (bytes == std::numeric_limits<std::size_t>::max()) {
            *TraceOptions::instance().out << std::string(4*nest,' ') << "[" << name_ << " (dealloc_async)] { pointer:" << p << ", stream:" << stream.value() << " }" << std::endl;
        }
        else {
            *TraceOptions::instance().out << std::string(4*nest,' ') << "[" << name_ << " (dealloc_async)] { pointer:" << p << ", bytes:" << bytes_to_string(bytes) << ", alignment:" << alignment << ", stream:" << stream.value() << " }" << std::endl;
        }
    }
    if (async_mr) {
        async_mr->deallocate_async(p, bytes, alignment, stream);
    }
    else {
        mr_->deallocate(p, bytes, alignment);
    }
    nest--;
}



}
