/*
 * (C) Copyright 2024- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "DeviceMemoryResource.h"

#include <iostream>
#include <memory>

#include "hic/hic.h"

#include "pluto/pluto_config.h"
#include "pluto/runtime.h"
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

constant_init<DeviceMemoryResource>  device_res{};

async_memory_resource* device_resource() {
    return &device_res.obj;
}

memory_pool_resource* device_pool_resource() {
    static MemoryPoolResource resource(device_resource());
    return &resource;
}

void* DeviceMemoryResource::do_allocate(std::size_t bytes, alignment_t) {
    void* ptr;
    if constexpr (PLUTO_HAVE_HIC) {
        if (devices()) {
            HIC_CALL(hicMalloc(&ptr, bytes));
        }
        else {
            ptr = new_delete_resource()->allocate(bytes, alignment);
        }
    }
    else {
        ptr = new_delete_resource()->allocate(bytes, alignment);
    }
    if constexpr (LOG) {
        std::cout << "               + hicMalloc(ptr:" << ptr << ", bytes:" << bytes << ");" << std::endl;
    }

    if (trace::enabled()) {
        std::string_view label = pluto::get_label();
        if (label.empty()) {
            trace::out << "PLUTO_TRACE " << "pluto::device_resource::allocate(bytes="<<trace::format_bytes(bytes)<<", alignment="<<alignment<<") -> ptr=" << ptr << '\n';
        }
        else {
            trace::out << "PLUTO_TRACE " << "pluto::device_resource::allocate(label="<<label<<", bytes="<<trace::format_bytes(bytes)<<", alignment="<<alignment<<") -> ptr=" << ptr << '\n';
        }
    }

    return ptr;
}

void DeviceMemoryResource::do_deallocate(void* ptr, std::size_t bytes, alignment_t) {

    if (trace::enabled()) {
        std::string_view label = pluto::get_label();
        if (label.empty()) {
            trace::out << "PLUTO_TRACE " << "pluto::device_resource::deallocate(ptr="<<ptr<<", bytes="<<trace::format_bytes(bytes)<<", alignment="<<alignment<<")" << std::endl;
        }
        else {
            trace::out << "PLUTO_TRACE " << "pluto::device_resource::deallocate(label="<<label<<", ptr="<<ptr<<", bytes="<<trace::format_bytes(bytes)<<", alignment="<<alignment<<")" << std::endl;
        }
    }

    if constexpr (PLUTO_HAVE_HIC) {
        if (devices()) {
            HIC_CALL(hicFree(ptr));
        }
        else {
            new_delete_resource()->deallocate(ptr, bytes, alignment);
        }
    }
    else {
        new_delete_resource()->deallocate(ptr, bytes, alignment);
    }
    if constexpr (LOG) {
        std::cout << "               - hicFree(ptr:" << ptr << ", bytes:" << bytes << ");" << std::endl;
    }
}


void* DeviceMemoryResource::do_allocate_async(std::size_t bytes, alignment_t, stream_view s) {
    void* ptr;
    if constexpr (PLUTO_HAVE_HIC) {
        HIC_CALL(hicMallocAsync(&ptr, bytes, s.value<hicStream_t>()));
    }
    else {
        ptr = new_delete_resource()->allocate(bytes, alignment);
    }
    if constexpr (LOG) {
        std::cout << "               + hicMallocAsync(ptr:" << ptr << ", bytes:" << bytes << ", stream:" << s.value()
                  << ");" << std::endl;
    }

    if (trace::enabled()) {
        std::string_view label = pluto::get_label();
        if (label.empty()) {
            trace::out << "PLUTO_TRACE " << "pluto::device_resource::allocate_async(bytes="<<trace::format_bytes(bytes)<<", alignment="<<alignment<<", stream="<< s.value()<<") -> ptr=" << ptr << '\n';
        }
        else {
            trace::out << "PLUTO_TRACE " << "pluto::device_resource::allocate_async(label="<<label<<", bytes="<<trace::format_bytes(bytes)<<", alignment="<<alignment<<", stream="<< s.value()<<") -> ptr=" << ptr << '\n';
        }
    }

    return ptr;
}

void DeviceMemoryResource::do_deallocate_async(void* ptr, std::size_t bytes, alignment_t, stream_view s) {
    if (trace::enabled()) {
        std::string_view label = pluto::get_label();
        if (label.empty()) {
            trace::out << "PLUTO_TRACE " << "pluto::device_resource::deallocate_async(ptr="<<ptr<<", bytes="<<trace::format_bytes(bytes)<<", alignment="<<alignment<<", stream="<<s.value()<<")" << std::endl;
        }
        else {
            trace::out << "PLUTO_TRACE " << "pluto::device_resource::deallocate_async(label="<<label<<", ptr="<<ptr<<", bytes="<<trace::format_bytes(bytes)<<", alignment="<<alignment<<", stream="<<s.value()<<")" << std::endl;
        }
    }

    if constexpr (PLUTO_HAVE_HIC) {
        HIC_CALL(hicFreeAsync(ptr, s.value<hicStream_t>()));
    }
    else {
        new_delete_resource()->deallocate(ptr, bytes, alignment);
    }
    if constexpr (LOG) {
        std::cout << "               - hicFreeAsync(ptr:" << ptr << ", bytes:" << bytes << ", stream:" << s.value()
                  << ");" << std::endl;
    }
}

bool DeviceMemoryResource::do_is_equal(const memory_resource& other) const noexcept {
    if (this == &other) {
        return true;
    }
    return false;
}

// ---------------------------------------------------------------------------------------------------------

}  // namespace pluto
