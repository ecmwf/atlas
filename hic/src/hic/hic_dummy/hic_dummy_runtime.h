/*
 * (C) Copyright 2024- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "hic/hic_dummy/dummyShouldNotBeCalled.h"

#define DUMMY_SHOULD_NOT_BE_CALLED(SYMBOL) dummyShouldNotBeCalled(#SYMBOL)
#define DUMMY_FUNCTION(SYMBOL)                     \
    template <typename... Args>                    \
    inline dummyError_t dummy##SYMBOL(Args&&...) { \
        DUMMY_SHOULD_NOT_BE_CALLED(hic##SYMBOL);   \
        return dummyError_t{0};                    \
    }
#define DUMMY_VALUE(SYMBOL) constexpr int dummy##SYMBOL = 0;

namespace {

using dummyError_t  = int;
using dummyEvent_t  = void*;
using dummyHostFn_t = void*;
using dummyStream_t = void*;

inline const char* dummyGetErrorString(dummyError_t) {
    DUMMY_SHOULD_NOT_BE_CALLED(hicGetErrorString);
}

inline dummyError_t dummyGetLastError(void) {
    DUMMY_SHOULD_NOT_BE_CALLED(hicGetLastError);
}

inline dummyError_t dummyPeekAtLastError(void) {
    DUMMY_SHOULD_NOT_BE_CALLED(hicPeekAtLastError);
}

struct dummyPointerAttributes {
    int type{0};
    int device{-2};
    void* hostPointer{nullptr};
    void* devicePointer{nullptr};
};

DUMMY_FUNCTION(DeviceSynchronize)
DUMMY_FUNCTION(Free)
DUMMY_FUNCTION(FreeAsync)
DUMMY_FUNCTION(GetDeviceCount)
DUMMY_FUNCTION(GetErrorString)
DUMMY_FUNCTION(GetLastError)
DUMMY_FUNCTION(HostGetDevicePointer)
DUMMY_FUNCTION(HostRegister)
DUMMY_FUNCTION(HostUnregister)
DUMMY_FUNCTION(Malloc)
DUMMY_FUNCTION(MallocAsync)
DUMMY_FUNCTION(MallocManaged)
DUMMY_FUNCTION(Memcpy)
DUMMY_FUNCTION(Memcpy2D)
DUMMY_FUNCTION(MemcpyAsync)
DUMMY_FUNCTION(Memcpy2DAsync)
DUMMY_FUNCTION(MemPrefetchAsync)
DUMMY_FUNCTION(StreamCreate)
DUMMY_FUNCTION(StreamDestroy)
DUMMY_FUNCTION(StreamSynchronize)
DUMMY_FUNCTION(PointerGetAttributes)

DUMMY_VALUE(CpuDeviceId)
DUMMY_VALUE(HostRegisterMapped)
DUMMY_VALUE(MemoryTypeDevice)
DUMMY_VALUE(MemoryTypeHost)
DUMMY_VALUE(MemoryTypeUnregistered)
DUMMY_VALUE(MemoryTypeManaged)
DUMMY_VALUE(MemcpyDeviceToHost)
DUMMY_VALUE(MemcpyHostToDevice)
DUMMY_VALUE(Success)

}  // namespace

#undef DUMMY_FUNCTION
#undef DUMMY_VALUE
#undef DUMMY_SHOULD_NOT_BE_CALLED