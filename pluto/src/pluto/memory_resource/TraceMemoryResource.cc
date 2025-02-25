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

#include <iomanip>
#include <limits>
#include <sstream>

#include "pluto/stream.h"

namespace pluto {

static std::string bytes_to_string(std::size_t bytes) {
    constexpr double KB = 1024.;
    constexpr double MB = 1024. * KB;
    constexpr double GB = 1024. * MB;
    double b            = static_cast<double>(bytes);

    std::stringstream ss;
    ss << std::setprecision(2) << std::fixed;
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
    return ss.str();
}


int TraceMemoryResource::nest = 0;

void* TraceMemoryResource::do_allocate(std::size_t bytes, std::size_t alignment) {
    nest++;
    if (trace_enabled()) {
        auto label = get_label();
        trace() << std::string(4 * nest, ' ') << "[" << name_ << " (alloc)] { ";
        if (not label.empty()) {
            trace() << "label:" << label << ", ";
        }
        trace() << "bytes:" << bytes_to_string(bytes)
                << ", alignment:" << alignment
                << " }" << std::endl;
    }
    auto ptr = upstream_resource()->allocate(bytes, alignment);
    nest--;
    return ptr;
}

void TraceMemoryResource::do_deallocate(void* p, std::size_t bytes, std::size_t alignment) {
    nest++;
    if (trace_enabled()) {
        auto label = get_label();
        trace() << std::string(4 * nest, ' ') << "[" << name_ << " (dealloc)] { ";
        if (not label.empty()) {
            trace() << "label:" << label << ", ";
        }
        trace() << "pointer:" << p;
        if (bytes != std::numeric_limits<std::size_t>::max()) {
            trace() << ", bytes:" << bytes_to_string(bytes)
                    << ", alignment:" << alignment;
        }
        trace() << " }" << std::endl;
    }
    upstream_resource()->deallocate(p, bytes, alignment);
    nest--;
}

void* TraceMemoryResource::do_allocate_async(std::size_t bytes, std::size_t alignment, stream_view s) {
    nest++;
    auto* async_mr = dynamic_cast<async_memory_resource*>(upstream_resource());

    if (trace_enabled()) {
        auto label = get_label();
        trace() << std::string(4 * nest, ' ') << "[" << name_ << " (alloc_async)] { ";
        if (not label.empty()) {
            trace() << "label:" << label << ", ";
        }
        trace() << "bytes:" << bytes_to_string(bytes)
                << ", alignment:" << alignment
                << ", stream:" << s.value() 
                << " }" << std::endl;
    }
    void* ptr;
    if (async_mr) {
        ptr = async_mr->allocate_async(bytes, alignment, s);
    }
    else {
        ptr = upstream_resource()->allocate(bytes, alignment);
    }
    nest--;
    return ptr;
}

void TraceMemoryResource::do_deallocate_async(void* p, std::size_t bytes, std::size_t alignment, stream_view s) {
    nest++;
    auto* async_mr = dynamic_cast<async_memory_resource*>(upstream_resource());
    if (trace_enabled()) {
        auto label = get_label();
        if (bytes == std::numeric_limits<std::size_t>::max()) {
            trace() << std::string(4 * nest, ' ') << "[" << name_ << " (dealloc_async)] { pointer:" << p
                    << ", stream:" << s.value() << " }" << std::endl;
        }
        else {
            trace() << std::string(4 * nest, ' ') << "[" << name_ << " (dealloc_async)] { pointer:" << p
                    << ", bytes:" << bytes_to_string(bytes) << ", alignment:" << alignment << ", stream:" << s.value()
                    << " }" << std::endl;
        }
    }
    if (async_mr) {
        async_mr->deallocate_async(p, bytes, alignment, s);
    }
    else {
        upstream_resource()->deallocate(p, bytes, alignment);
    }
    nest--;
}


}  // namespace pluto
