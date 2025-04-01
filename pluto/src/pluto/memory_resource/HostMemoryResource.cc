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

#include "pluto/alignment.h"
#include "pluto/pluto_config.h"
#include "pluto/trace.h"
#include "pluto/memory.h"

#include "MemoryPoolResource.h"

#define LOG PLUTO_DEBUGGING

namespace pluto {

class HostMemoryResource : public memory_resource {
public:
    HostMemoryResource() = default;

    ~HostMemoryResource() = default;

private:
    void* do_allocate(std::size_t bytes, std::size_t alignment) override;
    void do_deallocate(void* ptr, std::size_t bytes, std::size_t alignment) override;

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

constant_init<HostMemoryResource>  host_res{};

memory_resource* host_resource() {
    // Never destroyed due to constant_init!
    return &host_res.obj;
}

memory_pool_resource* host_pool_resource() {
    // Never destroyed due to constant_init!
    static constant_init<MemoryPoolResource> host_pool_res{host_resource(), "pluto::host_pool_resource", &pluto::memory::host_pool};
    return &host_pool_res.obj;
}


void* HostMemoryResource::do_allocate(std::size_t bytes, std::size_t alignment) {
    alignment = std::max(alignment, default_alignment());
    auto* ptr = new_delete_resource()->allocate(bytes, alignment);

    memory::host.allocate(bytes);
    if (trace::enabled()) {
        trace::log::allocate(get_label(), ptr, bytes, alignment, "pluto::host_resource", &memory::host);
    }

    return ptr;
}

void HostMemoryResource::do_deallocate(void* ptr, std::size_t bytes, std::size_t alignment) {
    alignment = std::max(alignment, default_alignment());

    memory::host.deallocate(bytes);
    if (trace::enabled()) {
        trace::log::deallocate(get_label(), ptr, bytes, alignment, "pluto::host_resource", &memory::host);
    }

    new_delete_resource()->deallocate(ptr, bytes, alignment);
}

bool HostMemoryResource::do_is_equal(const memory_resource& other) const noexcept {
    return (this == &other);
}

// ---------------------------------------------------------------------------------------------------------

}  // namespace pluto
