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

#if ATLAS_HAVE_CUDA
#include "hic/hic.h"
#endif

namespace atlas {
namespace util {

//------------------------------------------------------------------------------
namespace detail {
//------------------------------------------------------------------------------

void allocate_cudamanaged(void** ptr, size_t size) {
#if ATLAS_HAVE_CUDA
    hicError_t err = hicMallocManaged(ptr, size);
    if (err != hicSuccess)
        throw_AssertionFailed("failed to allocate GPU memory", Here());
#else
    *ptr = malloc(size);
#endif
}

void deallocate_cudamanaged(void* ptr) {
#if ATLAS_HAVE_CUDA
    hicError_t err = hicDeviceSynchronize();
    if (err != hicSuccess)
        throw_AssertionFailed("failed to synchronize memory", Here());

    err = hicFree(ptr);
    // The following throws an invalid device memory
    if (err != hicSuccess)
        throw_AssertionFailed("failed to free GPU memory", Here());
#else
    free(ptr);
#endif
}

void allocate_cuda(void** ptr, size_t size) {
#if ATLAS_HAVE_CUDA
    hicError_t err = hicMalloc(ptr, size);
    if (err != hicSuccess)
        throw_AssertionFailed("failed to allocate GPU memory", Here());
#else
    *ptr = malloc(size);
#endif
}

void deallocate_cuda(void* ptr) {
    deallocate_cudamanaged(ptr);
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
