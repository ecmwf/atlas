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
#include <string_view>

#include "hic/hic.h"

#include "pluto/pluto_config.h"
#include "pluto/runtime.h"
#include "pluto/memory.h"
#include "pluto/trace.h"

#include "MemoryPoolResource.h"

#define LOG PLUTO_DEBUGGING

namespace pluto {

// --------------------------------------------------------------------------------------------------------

class PinnedMemoryResource : public memory_resource {
public:
    PinnedMemoryResource() = default;

    void pin(void* ptr, std::size_t bytes);
    void unpin(void* ptr, std::size_t bytes);

    void* do_allocate(std::size_t bytes, std::size_t alignment) override;

    void do_deallocate(void* ptr, std::size_t bytes, std::size_t alignment) override;

    bool do_is_equal(const memory_resource& other) const noexcept override;
};

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

    pin(ptr, bytes);

    memory::pinned.allocate(bytes);
    if (trace::enabled()) {
        trace::log::allocate(get_label(), ptr, bytes, alignment, "pluto::pinned_resource", &memory::pinned);
    }
    return ptr;
}

void PinnedMemoryResource::do_deallocate(void* ptr, std::size_t bytes, std::size_t alignment) {
    alignment = std::max(alignment, default_alignment());

    memory::pinned.deallocate(bytes);
    if (trace::enabled()) {
        trace::log::deallocate(get_label(), ptr, bytes, alignment, "pluto::pinned_resource", &memory::pinned);
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

    template<typename ...Args>
    explicit constexpr constant_init(Args... args) : obj(args...) { }

    ~constant_init() { /* do nothing, union member is not destroyed */ }
};
}

constant_init<PinnedMemoryResource>  pinned_res{};

memory_resource* pinned_resource() {
    // Never destroyed due to constant_init!
    return &pinned_res.obj;
}

memory_pool_resource* pinned_pool_resource() {
    // Never destroyed due to constant_init!
    static constant_init<MemoryPoolResource> pinned_pool_res{pinned_resource(), "pluto::pinned_pool_resource", &pluto::memory::pinned_pool};
    return &pinned_pool_res.obj;
}

// --------------------------------------------------------------------------------------------------------

}  // namespace pluto
