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

#include "pluto/alignment.h"
#include "pluto/pluto_config.h"
#include "pluto/runtime.h"
#include "pluto/stream.h"
#include "pluto/memory.h"
#include "pluto/trace.h"

#include "MemoryPoolResource.h"

#define LOG PLUTO_DEBUGGING

namespace pluto {

class DeviceMemoryResource : public async_memory_resource {
public:
    using alignment_t                      = std::size_t;
    static constexpr alignment_t alignment = default_alignment();

    DeviceMemoryResource() = default;

    ~DeviceMemoryResource() = default;

private:
    void* do_allocate(std::size_t bytes, alignment_t) override;
    void do_deallocate(void* ptr, std::size_t bytes, std::size_t alignment) override;

    bool do_is_equal(const memory_resource& other) const noexcept override;

    void* do_allocate_async(std::size_t bytes, alignment_t, stream_view) override;
    void do_deallocate_async(void* ptr, std::size_t bytes, std::size_t alignment, stream_view) override;
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

constant_init<DeviceMemoryResource>  device_res{};

async_memory_resource* device_resource() {
    // Never destroyed due to constant_init!
    return &device_res.obj;
}

memory_pool_resource* device_pool_resource() {
    // Never destroyed due to constant_init!
    static constant_init<MemoryPoolResource> device_pool_res{device_resource(), "pluto::device_pool_resource", &pluto::memory::device_pool};
    return &device_pool_res.obj;
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

    memory::device.allocate(bytes);
    if (trace::enabled()) {
        trace::log::allocate(get_label(), ptr, bytes, alignment, "pluto::device_resource", &memory::device);
    }

    return ptr;
}

void DeviceMemoryResource::do_deallocate(void* ptr, std::size_t bytes, alignment_t) {

    memory::device.deallocate(bytes);
    if (trace::enabled()) {
        trace::log::deallocate(get_label(), ptr, bytes, alignment, "pluto::device_resource", &memory::device);
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

    memory::device.allocate(bytes);
    if (trace::enabled()) {
        trace::log::allocate_async(get_label(), ptr, bytes, alignment, s.value(), "pluto::device_resource", &memory::device);
    }

    return ptr;
}

void DeviceMemoryResource::do_deallocate_async(void* ptr, std::size_t bytes, alignment_t, stream_view s) {
    memory::device.deallocate(bytes);
    if (trace::enabled()) {
        trace::log::deallocate_async(get_label(), ptr, bytes, alignment, s.value(), "pluto::device_resource", &memory::device);
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
    return (this == &other);
}

// ---------------------------------------------------------------------------------------------------------

}  // namespace pluto
