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
#include "pluto/stream.h"

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
    void* ptr = upstream_->allocate(bytes, alignment);
    pin(ptr, bytes);
    return ptr;
}

void PinnedMemoryResource::do_deallocate(void* ptr, std::size_t bytes, std::size_t alignment) {
    unpin(ptr, bytes);
    upstream_->deallocate(ptr, bytes, alignment);
}

bool PinnedMemoryResource::do_is_equal(const memory_resource& other) const noexcept {
    return upstream_->is_equal(other);
}


// --------------------------------------------------------------------------------------------------------

memory_resource* pinned_resource() {
    static PinnedMemoryResource resource;
    return &resource;
}

memory_pool_resource* pinned_pool_resource() {
    static MemoryPoolResource resource(pinned_resource());
    return &resource;
}

}  // namespace pluto
