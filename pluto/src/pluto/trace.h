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

#include <cstddef>
#include <iostream>
#include <string>
#include <string_view>

namespace pluto {
    class memory_tracker;
}

namespace pluto::trace {

struct Options {
    bool enabled{false};
    std::ostream* out{&std::cout};

    Options() = default;
    Options(bool _enabled) : enabled{_enabled} {}
    Options(const Options& other) : 
        enabled{other.enabled},
        out{other.out} {}
    Options& operator=(Options& other) {
        enabled = other.enabled;
        out = other.out;
        return *this;
    }
};

Options& options();

struct OutputStream {
    template <typename T>
    std::ostream& operator<<(const T& value) {
        std::ostream& out = *options().out;
        out << value;
        return out;
    }
};
extern OutputStream out;


inline bool enabled() {
    return options().enabled;
}
inline bool enable(bool value = true) {
    auto& opts = options();
    bool previously_enabled = opts.enabled;
    opts.enabled = value;
    return previously_enabled;
}
inline void set(std::ostream& out) {
    options().out = &out;
}
inline void set_options(const Options& opts) {
    options().enabled = opts.enabled;
    options().out = opts.out;
}

std::string format_bytes(std::size_t bytes);

namespace log {
void allocate(std::string_view label, void* ptr, std::size_t bytes, std::size_t alignment, std::string_view resource_name, memory_tracker* memory_tracker);
void deallocate(std::string_view label, void* ptr, std::size_t bytes, std::size_t alignment, std::string_view resource_name, memory_tracker* memory_tracker);
void allocate_async(std::string_view label, void* ptr, std::size_t bytes, std::size_t alignment, void* stream, std::string_view resource_name, memory_tracker* memory_tracker);
void deallocate_async(std::string_view label, void* ptr, std::size_t bytes, std::size_t alignment, void* stream, std::string_view resource_name, memory_tracker* memory_tracker);
}

}  // namespace pluto::trace
