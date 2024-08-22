/*
 * (C) Copyright 2024- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "PinnedMemoryResource.h"

#include <iostream>
#include <mutex>

#include "hic/hic.h"

#include "pluto/pluto_config.h"
#include "pluto/runtime.h"
#include "pluto/trace.h"

#include "MemoryPoolResource.h"

#define LOG PLUTO_DEBUGGING

namespace pluto {

// --------------------------------------------------------------------------------------------------------

void PinnedMemoryResource::pin(void* ptr, std::size_t bytes) {
    if (LOG) {
        std::cout << "         + PIN: hicHostRegister(ptr:" << ptr << ", bytes:" << bytes << ");"
                  << std::endl;  //, hicHostRegisterMapped);" << std::endl;
    }
    if constexpr (PLUTO_HAVE_HIC) {
        if (devices()) {
            HIC_CALL(hicHostRegister(ptr, bytes, hicHostRegisterMapped));
        }
    }
}

// --------------------------------------------------------------------------------------------------------

void PinnedMemoryResource::unpin(void* ptr, std::size_t bytes) {
    if (LOG) {
        std::cout << "         - UNPIN: hicHostUnRegister(ptr:" << ptr << ", bytes:" << bytes << ");" << std::endl;
    }
    if constexpr (PLUTO_HAVE_HIC) {
        if (devices()) {
            HIC_CALL(hicHostUnregister(ptr));
        }
    }
}

void* PinnedMemoryResource::do_allocate(std::size_t bytes, std::size_t alignment) {
    alignment = std::max(alignment, default_alignment());
    auto* ptr = new_delete_resource()->allocate(bytes, alignment);

    if (trace::enabled()) {
        std::string_view label = pluto::get_label();
        if (label.empty()) {
            trace::out << "PLUTO_TRACE " << "pluto::pinned_resource::allocate(bytes="<<trace::format_bytes(bytes)<<", alignment="<<alignment<<") -> ptr=" << ptr << '\n';
        }
        else {
            trace::out << "PLUTO_TRACE " << "pluto::pinned_resource::allocate(label="<<label<<", bytes="<<trace::format_bytes(bytes)<<", alignment="<<alignment<<") -> ptr=" << ptr << '\n';
        }
    }

    pin(ptr, bytes);
    return ptr;
}

void PinnedMemoryResource::do_deallocate(void* ptr, std::size_t bytes, std::size_t alignment) {
    alignment = std::max(alignment, default_alignment());

    if (trace::enabled()) {
        std::string_view label = pluto::get_label();
        if (label.empty()) {
            trace::out << "PLUTO_TRACE " << "pluto::pinned_resource::deallocate(ptr="<<ptr<<", bytes="<<trace::format_bytes(bytes)<<", alignment="<<alignment<<")" << std::endl;
        }
        else {
            trace::out << "PLUTO_TRACE " << "pluto::pinned_resource::deallocate(label="<<label<<", ptr="<<ptr<<", bytes="<<trace::format_bytes(bytes)<<", alignment="<<alignment<<")" << std::endl;
        }
    }

    unpin(ptr, bytes);
    new_delete_resource()->deallocate(ptr, bytes, alignment);
}

bool PinnedMemoryResource::do_is_equal(const memory_resource& other) const noexcept {
    return (this == &other);
}

// --------------------------------------------------------------------------------------------------------

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

constant_init<PinnedMemoryResource>  pinned_res{};

memory_resource* pinned_resource() {
    return &pinned_res.obj;
}

memory_pool_resource* pinned_pool_resource() {
    static MemoryPoolResource resource(pinned_resource());
    return &resource;
}

}  // namespace pluto
