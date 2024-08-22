/*
 * (C) Copyright 2024- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "HostMemoryResource.h"

#include <iostream>
#include <memory>

#include "pluto/pluto_config.h"
#include "pluto/stream.h"
#include "pluto/trace.h"

#include "MemoryPoolResource.h"

#define LOG PLUTO_DEBUGGING

namespace pluto {

namespace {
template<typename T>
struct constant_init {
    union {
        T obj;
    };
    constexpr constant_init() : obj() { }

    template<typename U>
    explicit constexpr constant_init(U arg) : obj(arg) { }

    ~constant_init() { /* do nothing, union member is not destroyed */ }
};
}

constant_init<HostMemoryResource>  host_res{};

async_memory_resource* host_resource() {
    return &host_res.obj;
}

memory_pool_resource* host_pool_resource() {
    static MemoryPoolResource resource(host_resource());
    return &resource;
}


void* HostMemoryResource::do_allocate(std::size_t bytes, std::size_t alignment) {
    alignment = std::max(alignment, default_alignment());
    auto* ptr = new_delete_resource()->allocate(bytes, alignment);

    if (trace::enabled()) {
        std::string_view label = pluto::get_label();
        if (label.empty()) {
            trace::out << "PLUTO_TRACE " << "pluto::host_resource::allocate(bytes="<<trace::format_bytes(bytes)<<", alignment="<<alignment<<") -> ptr=" << ptr << '\n';
        }
        else {
            trace::out << "PLUTO_TRACE " << "pluto::host_resource::allocate(label="<<label<<", bytes="<<trace::format_bytes(bytes)<<", alignment="<<alignment<<") -> ptr=" << ptr << '\n';
        }
    }

    return ptr;
}

void HostMemoryResource::do_deallocate(void* ptr, std::size_t bytes, std::size_t alignment) {
    alignment = std::max(alignment, default_alignment());

    if (trace::enabled()) {
        std::string_view label = pluto::get_label();
        if (label.empty()) {
            trace::out << "PLUTO_TRACE " << "pluto::host_resource::deallocate(ptr="<<ptr<<", bytes="<<trace::format_bytes(bytes)<<", alignment="<<alignment<<")" << std::endl;
        }
        else {
            trace::out << "PLUTO_TRACE " << "pluto::host_resource::deallocate(label="<<label<<", ptr="<<ptr<<", bytes="<<trace::format_bytes(bytes)<<", alignment="<<alignment<<")" << std::endl;
        }
    }

    new_delete_resource()->deallocate(ptr, bytes, alignment);
}


void* HostMemoryResource::do_allocate_async(std::size_t bytes, std::size_t alignment, stream_view) {
    return new_delete_resource()->allocate(bytes, std::max(alignment, default_alignment()));
}

void HostMemoryResource::do_deallocate_async(void* ptr, std::size_t bytes, std::size_t alignment, stream_view) {
    new_delete_resource()->deallocate(ptr, bytes, std::max(alignment, default_alignment()));
}

bool HostMemoryResource::do_is_equal(const memory_resource& other) const noexcept {
    return (this == &other);
}

// ---------------------------------------------------------------------------------------------------------

}  // namespace pluto
