/*
* (C) Copyright 2013 ECMWF.
*
* This software is licensed under the terms of the Apache Licence Version 2.0
* which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
* In applying this licence, ECMWF does not waive the privileges and immunities
* granted to it by virtue of its status as an intergovernmental organisation nor
* does it submit to any jurisdiction.
*/


#include "Allocate.h"

#include "eckit/log/CodeLocation.h"

#include "atlas/library/config.h"
#include "atlas/runtime/Exception.h"

#include "hic/hic.h"

namespace atlas {
namespace util {

//------------------------------------------------------------------------------
namespace detail {
//------------------------------------------------------------------------------

void allocate_managed(void** ptr, size_t size) {
    if constexpr (not ATLAS_HAVE_GPU) {
        return allocate_host(ptr, size);
    }
    HIC_CALL(hicMallocManaged(ptr, size));
}

void deallocate_managed(void* ptr) {
    if constexpr (not ATLAS_HAVE_GPU) {
        return deallocate_host(ptr);
    }
    HIC_CALL(hicDeviceSynchronize());
    HIC_CALL(hicFree(ptr));
}

void allocate_device(void** ptr, size_t size) {
    if constexpr (not ATLAS_HAVE_GPU) {
        return allocate_host(ptr, size);
    }
    HIC_CALL(hicMalloc(ptr, size));
}

void deallocate_device(void* ptr) {
    if constexpr (not ATLAS_HAVE_GPU) {
        return deallocate_host(ptr);
    }
    HIC_CALL(hicDeviceSynchronize());
    HIC_CALL(hicFree(ptr));
}

void allocate_host(void** ptr, size_t size) {
    *ptr = malloc(size);
}

void deallocate_host(void* ptr) {
    free(ptr);
}

//------------------------------------------------------------------------------
}  // namespace detail
//------------------------------------------------------------------------------

extern "C" {
void atlas__allocate_managedmem_double(double*& a, size_t N) {
    allocate_managedmem(a, N);
}
void atlas__allocate_managedmem_float(float*& a, size_t N) {
    allocate_managedmem(a, N);
}
void atlas__allocate_managedmem_int(int*& a, size_t N) {
    allocate_managedmem(a, N);
}
void atlas__allocate_managedmem_long(long*& a, size_t N) {
    allocate_managedmem(a, N);
}
void atlas__deallocate_managedmem(void*& a) {
    delete_managedmem(a);
}
}

//------------------------------------------------------------------------------

}  // namespace util
}  // namespace atlas
