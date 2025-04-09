/*
 * (C) Copyright 2024- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */
#pragma once

#include "hic/hic_namespace_macro.h"

#if HIC_BACKEND_CUDA
  #define HIC_BACKEND cuda
  #include <cuda_runtime.h>
#elif HIC_BACKEND_HIP
  #define HIC_BACKEND hip

  #pragma push_macro("DEPRECATED")
  #undef DEPRECATED
  #include <hip/hip_runtime.h>
  #pragma pop_macro("DEPRECATED")

#if HIP_VERSION_MAJOR < 6
enum hicMemoryType
{
    hicMemoryTypeUnregistered = 0,
    hicMemoryTypeHost         = 1,
    hicMemoryTypeDevice       = 2,
    hicMemoryTypeManaged      = 3
};
struct hicPointerAttributes {
    hicMemoryType type;
    int device;
    void* devicePointer;
    void* hostPointer;
};

[[nodiscard]] inline hipError_t hicPointerGetAttributes(hicPointerAttributes* attributes, const void* ptr) {
    hipPointerAttribute_t attr_;
    auto err = hipPointerGetAttributes(&attr_, ptr);
    if (err != hipSuccess) {
        attributes->device        = 0;
        attributes->devicePointer = nullptr;
        attributes->hostPointer   = nullptr;
        return err;
    }
    const auto& type = attr_.memoryType;
    if (attr_.isManaged) {
        attributes->type = hicMemoryTypeManaged;
    }
    else if (type == hipMemoryTypeHost) {
        attributes->type = hicMemoryTypeHost;
    }
    else if (type == hipMemoryTypeDevice) {
        attributes->type = hicMemoryTypeDevice;
    }
    else if (type == hipMemoryTypeUnified) {
        attributes->type = hicMemoryTypeManaged;
    }
    else {
        attributes->type = hicMemoryTypeUnregistered;
    }
    attributes->device        = attr_.device;
    attributes->devicePointer = attr_.devicePointer;
    attributes->hostPointer   = attr_.hostPointer;
    return err;
};
#endif

#elif HIC_BACKEND_DUMMY
  #define HIC_BACKEND dummy
  #include "hic/hic_dummy/hic_dummy_runtime.h"
#else
  #error Unsupported hic backend. Please define HIC_BACKEND_CUDA or HIC_BACKEND_HIP or HIC_BACKEND_DUMMY
#endif

#define HIC_PREFIX hic
#define HIC_CONCAT_(A, B) A##B
#define HIC_CONCAT(A, B) HIC_CONCAT_(A, B)
#define HIC_SYMBOL(API) HIC_CONCAT(HIC_PREFIX, API)
#define HIC_BACKEND_SYMBOL(API) HIC_CONCAT(HIC_BACKEND, API)

#define HIC_TYPE(TYPE) using HIC_SYMBOL(TYPE) = HIC_BACKEND_SYMBOL(TYPE);

#define HIC_FUNCTION(FUNCTION)                                                   \
    template <typename... Args>                                                  \
    inline auto HIC_SYMBOL(FUNCTION)(Args && ... args)                           \
        -> decltype(HIC_BACKEND_SYMBOL(FUNCTION)(std::forward<Args>(args)...)) { \
        return HIC_BACKEND_SYMBOL(FUNCTION)(std::forward<Args>(args)...);        \
    }

#define HIC_VALUE(VALUE) constexpr decltype(HIC_BACKEND_SYMBOL(VALUE)) HIC_SYMBOL(VALUE) = HIC_BACKEND_SYMBOL(VALUE);

//------------------------------------------------
HIC_NAMESPACE_BEGIN
//------------------------------------------------

HIC_TYPE(HostFn_t)
HIC_TYPE(Error_t)
HIC_TYPE(Event_t)
HIC_TYPE(Stream_t)
#if !HIC_BACKEND_HIP
HIC_TYPE(PointerAttributes)
#elif HIP_VERSION_MAJOR >= 6
using HIC_SYMBOL(PointerAttributes) = HIC_BACKEND_SYMBOL(PointerAttribute_t);
#endif

HIC_FUNCTION(DeviceSynchronize)
HIC_FUNCTION(Free)
HIC_FUNCTION(FreeAsync)
HIC_FUNCTION(GetDeviceCount)
HIC_FUNCTION(GetErrorString)
HIC_FUNCTION(GetLastError)
HIC_FUNCTION(HostGetDevicePointer)
HIC_FUNCTION(LaunchHostFunc)
HIC_FUNCTION(HostRegister)
HIC_FUNCTION(HostUnregister)
HIC_FUNCTION(Malloc)
HIC_FUNCTION(MallocAsync)
HIC_FUNCTION(MallocManaged)
HIC_FUNCTION(Memcpy)
HIC_FUNCTION(Memcpy2D)
HIC_FUNCTION(MemcpyAsync)
HIC_FUNCTION(Memcpy2DAsync)
HIC_FUNCTION(MemPrefetchAsync)
HIC_FUNCTION(PeekAtLastError)
#if !HIC_BACKEND_HIP
HIC_FUNCTION(PointerGetAttributes)
#elif HIP_VERSION_MAJOR >= 6
HIC_FUNCTION(PointerGetAttributes)
#endif
HIC_FUNCTION(StreamCreate)
HIC_FUNCTION(StreamDestroy)
HIC_FUNCTION(StreamSynchronize)

HIC_VALUE(CpuDeviceId)
HIC_VALUE(HostRegisterMapped)
#if !HIC_BACKEND_HIP || (HIC_BACKEND_HIP && HIP_VERSION_MAJOR >= 6)
HIC_VALUE(MemoryTypeDevice)
HIC_VALUE(MemoryTypeHost)
HIC_VALUE(MemoryTypeUnregistered)
HIC_VALUE(MemoryTypeManaged)
#endif
HIC_VALUE(MemcpyDeviceToHost)
HIC_VALUE(MemcpyHostToDevice)
HIC_VALUE(Success)

#if HIC_BACKEND_CUDA
    constexpr hicError_t hicErrorDeinitialized = cudaErrorCudartUnloading;
#elif HIC_BACKEND_HIP
    constexpr hicError_t hicErrorDeinitialized = hipErrorDeinitialized;
#else
    constexpr hicError_t hicErrorDeinitialized = 4;
#endif

//------------------------------------------------
HIC_NAMESPACE_END
//------------------------------------------------

#undef HIC_FUNCTION
#undef HIC_TYPE
#undef HIC_VALUE
#undef HIC_CONCAT
#undef HIC_CONCAT_
#undef HIC_PREFIX
#undef HIC_SYMBOL
#undef HIC_BACKEND_SYMBOL
