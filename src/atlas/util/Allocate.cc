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

#include "pluto/pluto.h"

#include "atlas/library/config.h"
#include "atlas/runtime/Exception.h"

namespace atlas {
namespace util {

//------------------------------------------------------------------------------
namespace detail {
//------------------------------------------------------------------------------

void allocate_managed(void** ptr, size_t bytes) {
    if constexpr (not ATLAS_HAVE_GPU) {
        return allocate_host(ptr, bytes);
    }
    *ptr = pluto::managed_resource()->allocate(bytes, pluto::default_alignment());
}

void deallocate_managed(void* ptr, size_t bytes) {
    if constexpr (not ATLAS_HAVE_GPU) {
        return deallocate_host(ptr, bytes);
    }
    pluto::wait();
    pluto::managed_resource()->deallocate(ptr, bytes, pluto::default_alignment());
}

void allocate_device(void** ptr, size_t bytes) {
    if constexpr (not ATLAS_HAVE_GPU) {
        return allocate_host(ptr, bytes);
    }
    *ptr = pluto::device_resource()->allocate(bytes, pluto::default_alignment());
}

void deallocate_device(void* ptr, size_t bytes) {
    if constexpr (not ATLAS_HAVE_GPU) {
        return deallocate_host(ptr, bytes);
    }
    pluto::wait();
    pluto::device_resource()->deallocate(ptr, bytes, pluto::default_alignment());
}

void allocate_host(void** ptr, size_t bytes) {
    *ptr = malloc(bytes);
}

void deallocate_host(void* ptr, size_t /*bytes*/) {
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
void atlas__deallocate_managedmem_double(double*& a, size_t N) {
    delete_managedmem(a, N);
}
void atlas__deallocate_managedmem_float(float*& a, size_t N) {
    delete_managedmem(a, N);
}
void atlas__deallocate_managedmem_int(int*& a, size_t N) {
    delete_managedmem(a, N);
}
void atlas__deallocate_managedmem_long(long*& a, size_t N) {
    delete_managedmem(a, N);
}
}

//------------------------------------------------------------------------------

}  // namespace util
}  // namespace atlas
