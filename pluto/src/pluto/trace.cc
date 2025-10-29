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
#include <vector>

#include "pluto/memory.h"

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

namespace log {
void allocate(std::string_view label, const void* ptr, std::size_t bytes, std::size_t alignment, std::string_view resource_name, memory_tracker* memory_tracker) {
    out << "PLUTO_TRACE " << resource_name << "::allocate(";
    if (not label.empty()) {
        out << "label="<<label<<", ";
    }
    out << "bytes="<<format_bytes(bytes)<<", alignment="<<alignment<<") -> ptr=" << ptr;
    if (memory_tracker) {
        out << ", high_watermark="<<format_bytes(memory_tracker->high_watermark());
    }
    out << '\n';
}
void allocate_async(std::string_view label, const void* ptr, std::size_t bytes, std::size_t alignment, const void* stream, std::string_view resource_name, memory_tracker* memory_tracker) {
    out << "PLUTO_TRACE " << resource_name << "::allocate_async(";
    if (not label.empty()) {
        out << "label="<<label<<", ";
    }
    out << "bytes="<<format_bytes(bytes)<<", alignment="<<alignment<<", stream="<<stream<<") -> ptr=" << ptr;
    if (memory_tracker) {
        out << ", high_watermark="<<format_bytes(memory_tracker->high_watermark());
    }
    out << '\n';
}

void deallocate(std::string_view label, const void* ptr, std::size_t bytes, std::size_t alignment, std::string_view resource_name, memory_tracker*) {
    out << "PLUTO_TRACE " << resource_name << "::deallocate(";
    if (not label.empty()) {
        out << "label="<<label<<", ";
    }
    out << "ptr="<<ptr<<", bytes="<<format_bytes(bytes)<<", alignment="<<alignment<<")";
    out << '\n';
}
void deallocate_async(std::string_view label, const void* ptr, std::size_t bytes, std::size_t alignment, const void* stream, std::string_view resource_name, memory_tracker*) {
    out << "PLUTO_TRACE " << resource_name << "::deallocate_async(";
    if (not label.empty()) {
        out << "label="<<label<<", ";
    }
    out << "ptr="<<ptr<<", bytes="<<format_bytes(bytes)<<", alignment="<<alignment<<", stream="<<stream<<")";
    out << '\n';
}

void copy_host_to_device(std::string_view label, const void* dptr, const void* hptr, std::size_t bytes) {
    out << "PLUTO_TRACE copy_host_to_device(";
    if (not label.empty()) {
        out << "label="<<label<<", ";
    }
    out << "device_ptr="<<dptr<<", host_ptr="<<hptr<<", bytes="<<format_bytes(bytes)<<")";
    out << '\n';
}

void copy_device_to_host(std::string_view label, const void* hptr, const void* dptr, std::size_t bytes) {
    out << "PLUTO_TRACE copy_device_to_host(";
    if (not label.empty()) {
        out << "label="<<label<<", ";
    }
    out << "host_ptr="<<hptr<<", device_ptr="<<dptr<<", bytes="<<format_bytes(bytes)<<")";
    out << '\n';
}

void copy_host_to_host(std::string_view label, const void* dst_ptr, const void* src_ptr, std::size_t bytes) {
    out << "PLUTO_TRACE copy_host_to_host(";
    if (not label.empty()) {
        out << "label="<<label<<", ";
    }
    out << "dst_ptr="<<dst_ptr<<", src_ptr="<<src_ptr<<", bytes="<<format_bytes(bytes)<<")";
    out << '\n';
}

}


}  // namespace pluto::trace
