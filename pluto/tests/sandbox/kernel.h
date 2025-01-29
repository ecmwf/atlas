/*
 * (C) Copyright 2024- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <cstddef>
#include <stdio.h>

#include "hic/hic.h"

#include "pluto/pluto.h"

void set_on_device(double* x, std::size_t size, double value);

void set_on_device(const pluto::stream& stream, double* x, std::size_t size, double value);

#if HIC_COMPILER

template<typename F>
__global__ void applyKernel(std::size_t size, F f) {
   int idx = threadIdx.x + blockIdx.x * blockDim.x;
   if (idx < size) {
    // printf("idx %d / %d \n",idx,int(size));
    f(idx);
   }
}

template<typename F>
__global__ void applyKernel(F f) {
   int idx = threadIdx.x + blockIdx.x * blockDim.x;
   if (idx == 0) {
    f();
   }
}

template <typename F>
void launch_kernel( std::size_t size, F f) {
    printf("launch on device \n");
    constexpr std::size_t maxThreadsPerBlock = 256;
    int threadsPerBlock = std::min(size,maxThreadsPerBlock);
    int numBlocks = (size + threadsPerBlock - 1)/threadsPerBlock;
    // printf("threadsPerBlock = %d\n", threadsPerBlock);
    // printf("numBlocks = %d\n", numBlocks);
    applyKernel<<<numBlocks,threadsPerBlock>>>(size, f);
    HIC_CHECK_KERNEL_LAUNCH();
    HIC_CALL(hicDeviceSynchronize());
}

template <typename F>
void launch_kernel( F f) {
    printf("launch on device \n");
    applyKernel<<<1,1>>>(f);
    HIC_CHECK_KERNEL_LAUNCH();
    HIC_CALL(hicDeviceSynchronize());
}

template <typename F>
void launch_kernel(const pluto::stream& stream, std::size_t size, F f) {
    printf("launch on device async \n");
    constexpr std::size_t maxThreadsPerBlock = 256;
    int threadsPerBlock = std::min(size,maxThreadsPerBlock);
    int numBlocks = (size + threadsPerBlock - 1)/threadsPerBlock;
    // printf("threadsPerBlock = %d\n", threadsPerBlock);
    // printf("numBlocks = %d\n", numBlocks);
    applyKernel<<<numBlocks,threadsPerBlock,0,stream.value<hicStream_t>()>>>(size, f);
    HIC_CHECK_KERNEL_LAUNCH();
}

template <typename F>
void launch_kernel(const pluto::stream& stream, F f) {
    printf("launch on device async \n");
    applyKernel<<<1,1,0,stream.value<hicStream_t>()>>>(f);
    HIC_CHECK_KERNEL_LAUNCH();
}


// template <typename Derived, typename T>
// __global__ void apply_construct(T*& ptr) {
//     int idx = threadIdx.x + blockIdx.x * blockDim.x;
//     if (idx == 0) {
//         ptr = new Derived();
//     }
// }

// template <typename Derived, typename T>
// __device__ void apply_construct(T*& ptr) {
//     int idx = threadIdx.x + blockIdx.x * blockDim.x;
//     if (idx == 0) {
//         ptr = new Derived();
//     }
// }

template <typename Derived, typename T, typename... Args>
__global__ void apply_construct(T*& ptr, Args... args) {
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    if (idx == 0) {
        ptr = new Derived(std::forward<Args>(args)...);
    }
}

template <typename Derived, typename T, typename... Args>
void launch_construct( T*& ptr, Args... args) {
    printf("construct on device \n");
    apply_construct<Derived><<<1,1>>>(ptr, std::forward<Args>(args)...);
    HIC_CHECK_KERNEL_LAUNCH();
    HIC_CALL(hicDeviceSynchronize());
}


template <typename T>
__global__ void apply_destruct(T*& ptr) {
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    if (idx == 0) {
        delete ptr;
    }
}

template <typename T>
void launch_destruct( T*& ptr) {
    printf("destruct on device \n");
    apply_destruct<<<1,1>>>(ptr);
    HIC_CHECK_KERNEL_LAUNCH();
    HIC_CALL(hicDeviceSynchronize());
}

#else

template <typename F>
void launch_kernel(std::size_t size, F f) {
    printf("launch on host \n");
    for (size_t idx=0; idx<size; ++idx) {
        f(idx);
    }
}

template <typename F>
void launch_kernel(F f) {
    printf("launch on host \n");
    f();
}

template <typename F>
void launch_kernel(const pluto::stream& stream, std::size_t size, F f) {
    printf("launch on host \n");
    for (size_t idx=0; idx<size; ++idx) {
        f(idx);
    }
}

template <typename F>
void launch_kernel(const pluto::stream& stream, F f) {
    printf("launch on host \n");
    f();
}


#endif
