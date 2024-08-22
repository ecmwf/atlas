/*
 * (C) Copyright 2024- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "ManagedMemoryResource.h"

#include <iostream>

#include "hic/hic.h"
#include "pluto/pluto_config.h"
#include "pluto/runtime.h"
#include "pluto/stream.h"
#include "pluto/trace.h"

#include "MemoryPoolResource.h"
#include "TraceMemoryResource.h"

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

constant_init<ManagedMemoryResource>  managed_res{};

memory_resource* managed_resource() {
    return &managed_res.obj;
}

memory_pool_resource* managed_pool_resource() {
    static MemoryPoolResource resource(managed_resource());
    return &resource;
}

// ---------------------------------------------------------------------------------------------------------

void* ManagedMemoryResource::do_allocate(std::size_t bytes, alignment_t) {
    void* ptr = nullptr;
    if constexpr (PLUTO_HAVE_HIC) {
        if (devices()) {
            HIC_CALL(hicMallocManaged(&ptr, bytes));
        }
        else {
            ptr = new_delete_resource()->allocate(bytes, alignment);
        }
    }
    else {
        ptr = new_delete_resource()->allocate(bytes, alignment);
    }
    if constexpr (LOG) {
        std::cout << "               + hicMallocManaged(ptr:" << ptr << ", bytes:" << bytes << ")" << std::endl;
    }

    if (trace::enabled()) {
        std::string_view label = pluto::get_label();
        if (label.empty()) {
            trace::out << "PLUTO_TRACE " << "pluto::managed_resource::allocate(ptr="<<ptr<<", bytes="<<trace::format_bytes(bytes)<<", alignment="<<alignment<<")" << std::endl;
        }
        else {
            trace::out << "PLUTO_TRACE " << "pluto::managed_resource::allocate(label="<<label<<", ptr="<<ptr<<", bytes="<<trace::format_bytes(bytes)<<", alignment="<<alignment<<")" << std::endl;
        }
    }

    return ptr;
}

void ManagedMemoryResource::do_deallocate(void* ptr, std::size_t bytes, alignment_t) {
    if (trace::enabled()) {
        std::string_view label = pluto::get_label();
        if (label.empty()) {
            trace::out << "PLUTO_TRACE " << "pluto::managed_resource::deallocate(ptr="<<ptr<<", bytes="<<trace::format_bytes(bytes)<<", alignment="<<alignment<<")" << std::endl;
        }
        else {
            trace::out << "PLUTO_TRACE " << "pluto::managed_resource::deallocate(label="<<label<<", ptr="<<ptr<<", bytes="<<trace::format_bytes(bytes)<<", alignment="<<alignment<<")" << std::endl;
        }
    }

    if constexpr (LOG) {
        std::cout << "               - hicFree(ptr:" << ptr << ", bytes:" << bytes << ")" << std::endl;
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
}

bool ManagedMemoryResource::do_is_equal(const memory_resource& other) const noexcept {
    return (this == &other);
}

// ---------------------------------------------------------------------------------------------------------

}  // namespace pluto
