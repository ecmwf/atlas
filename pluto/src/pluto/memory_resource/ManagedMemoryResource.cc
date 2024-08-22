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

#include "pluto/pluto_config.h"
#include "pluto/offload/Stream.h"
#include "pluto/util/Runtime.h"
#include "hic/hic.h"

#include "MemoryPoolResource.h"
#include "TraceMemoryResource.h"

#define LOG PLUTO_DEBUGGING

namespace pluto {


memory_resource* managed_resource() {
    static ManagedMemoryResource resource;
    return &resource;
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
            HIC_CALL( hicMallocManaged(&ptr, bytes) );
        }
        else {
            ptr = std::pmr::new_delete_resource()->allocate(bytes, alignment);
        }
    }
    else {
        ptr = std::pmr::new_delete_resource()->allocate(bytes, alignment);
    }
    if constexpr (LOG) {
        std::cout << "               + hicMallocManaged(ptr:"<<ptr<<", bytes:"<< bytes <<")" << std::endl;
    }
    return ptr;
}

void ManagedMemoryResource::do_deallocate(void* ptr, std::size_t bytes, alignment_t) {
    if constexpr (LOG) {
        std::cout << "               - hicFree(ptr:"<<ptr<<", bytes:"<< bytes <<")" << std::endl;
    }
    if constexpr (PLUTO_HAVE_HIC) {
        if (devices()) {
            HIC_CALL( hicFree(ptr) );
        }
        else {
            std::pmr::new_delete_resource()->deallocate(ptr, bytes, alignment);
        }
    }
    else {
        std::pmr::new_delete_resource()->deallocate(ptr, bytes, alignment);
    }
}

bool ManagedMemoryResource::do_is_equal(const memory_resource& other) const noexcept {
    if (this == &other) { return true; }
    return false;
}

// ---------------------------------------------------------------------------------------------------------

}
