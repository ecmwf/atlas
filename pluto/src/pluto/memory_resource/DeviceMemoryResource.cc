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
#include "pluto/offload/Stream.h"

#include "MemoryPoolResource.h"

#define LOG PLUTO_DEBUGGING

namespace pluto {

memory_resource* device_resource() {
    static DeviceMemoryResource resource;
    return &resource;
}

memory_pool_resource* device_pool_resource() {
    static MemoryPoolResource resource(std::make_unique<DeviceMemoryResource>());
    return &resource;
}

void* DeviceMemoryResource::do_allocate(std::size_t bytes, alignment_t) {
    void* ptr;
    const auto& stream = get_default_stream();
    if (false) {//stream.value()) {
        if constexpr (PLUTO_HAVE_HIC) {
            HIC_CALL( hicMallocAsync(&ptr, bytes, stream.value<hicStream_t>() ) );
        }
        else {
            ptr = std::pmr::new_delete_resource()->allocate(bytes, alignment);
        }
        if constexpr (LOG) {
            std::cout << "               + hicMallocAsync(ptr:"<<ptr<<", bytes:"<< bytes <<", stream:" << stream.value() <<");" << std::endl;
        }
    }
    else {
        if constexpr (PLUTO_HAVE_HIC) {
            HIC_CALL( hicMalloc(&ptr, bytes) );
        }
        else {
            ptr = std::pmr::new_delete_resource()->allocate(bytes, alignment);
        }
        if constexpr (LOG) {
            std::cout << "               + hicMalloc(ptr:"<<ptr<<", bytes:"<< bytes <<");" << std::endl;
        }
    }
    return ptr;
}

void DeviceMemoryResource::do_deallocate(void* ptr, std::size_t bytes, alignment_t) {
    const auto& stream = get_default_stream();
    if (false) {//if (stream.value()) {
        if constexpr (PLUTO_HAVE_HIC) {
            HIC_CALL( hicFreeAsync(ptr, stream.value<hicStream_t>()) );
        }
        else {
            std::pmr::new_delete_resource()->deallocate(ptr, bytes, alignment);
        }
        if constexpr (LOG) {
            std::cout << "               - hicFreeAsync(ptr:"<<ptr<<", bytes:"<< bytes << ", stream:"<<stream.value()<<");" << std::endl;
        }
    }
    else {
        if constexpr (PLUTO_HAVE_HIC) {
            HIC_CALL( hicFree(ptr) );
        }
        else {
            std::pmr::new_delete_resource()->deallocate(ptr, bytes, alignment);
        }
        if constexpr (LOG) {
            std::cout << "               - hicFree(ptr:"<<ptr<<", bytes:"<< bytes << ");" << std::endl;
        }
    }
}

bool DeviceMemoryResource::do_is_equal(const memory_resource& other) const noexcept {
    if (this == &other) { return true; }
    return false;
}

// ---------------------------------------------------------------------------------------------------------

}
