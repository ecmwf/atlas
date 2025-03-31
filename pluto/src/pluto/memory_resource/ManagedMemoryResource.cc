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

#include "pluto/alignment.h"
#include "pluto/pluto_config.h"
#include "pluto/runtime.h"
#include "pluto/stream.h"
#include "pluto/memory.h"
#include "pluto/trace.h"

#include "MemoryPoolResource.h"

#define LOG PLUTO_DEBUGGING

namespace pluto {

class ManagedMemoryResource : public memory_resource {
public:
    using alignment_t                      = std::size_t;
    static constexpr alignment_t alignment = default_alignment();

    ManagedMemoryResource() = default;

protected:
    void* do_allocate(std::size_t bytes, alignment_t) override;
    void do_deallocate(void* ptr, std::size_t bytes, alignment_t) override;
    bool do_is_equal(const memory_resource& other) const noexcept override;
};

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

constant_init<ManagedMemoryResource>  managed_res{};

memory_resource* managed_resource() {
    // Never destroyed due to constant_init!
    return &managed_res.obj;
}

memory_pool_resource* managed_pool_resource() {
    // Never destroyed due to constant_init!
    static constant_init<MemoryPoolResource> managed_pool_res{managed_resource(), "pluto::managed_pool_resource", &pluto::memory::managed_pool};
    return &managed_pool_res.obj;
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

    memory::managed.allocate(bytes);
    if (trace::enabled()) {
        trace::log::allocate(get_label(), ptr, bytes, alignment, "pluto::managed_resource", &memory::managed);
    }

    return ptr;
}

void ManagedMemoryResource::do_deallocate(void* ptr, std::size_t bytes, alignment_t) {

    memory::managed.deallocate(bytes);
    if (trace::enabled()) {
        trace::log::deallocate(get_label(), ptr, bytes, alignment, "pluto::managed_resource", &memory::managed);
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
