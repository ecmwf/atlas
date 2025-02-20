#pragma once

#include <chrono>
#include <functional>
#include <iostream>

#include "hic/hic.h"
#include "pluto/pluto_config.h"

float run_timed(std::function<void()> function, int num_repeats = 1, int num_warmups = 0);

#if HIC_COMPILER
template <typename F>
__global__ void run_kernel_on_device(uint32_t n, F f) {
    const uint32_t idx{blockDim.x * blockIdx.x + threadIdx.x};
    const uint32_t stride{blockDim.x * gridDim.x};
    for (uint32_t i{idx}; i < n; i += stride) {
        f(i);
    }
}

template <typename F>
void launch_benchmark_kernel(uint32_t n, F f) {
    HIC_CALL(hicDeviceSynchronize());
    const uint32_t threads_per_block{1024};
    const uint32_t blocks_per_grid{32};
    run_kernel_on_device<<<blocks_per_grid, threads_per_block>>>(n, f);
    HIC_CALL(hicDeviceSynchronize());
}
#else
template <typename F>
void launch_benchmark_kernel(uint32_t n, F f) {
    for (uint32_t i{0}; i < n; ++i) {
        f(i);
    }
}
#endif


template <typename T>
inline void initialize_host_memory(T* h_buffer, uint32_t n, T value) {
    for (uint32_t i{0}; i < n; ++i) {
        h_buffer[i] = value;
    }
}

template <typename T>
inline bool verify_host_memory(T* h_buffer, uint32_t n, T value) {
    for (uint32_t i{0}; i < n; ++i) {
        if (h_buffer[i] != value) {
            return false;
        }
    }
    return true;
}

float launch_benchmark(float* d_output, float const* d_input_1, float const* d_input_2, uint32_t n);
float launch_benchmark(double* d_output, double const* d_input_1, double const* d_input_2, uint32_t n);

float touch_host_memory(float* h_buffer, uint32_t n);
float touch_host_memory(double* h_buffer, uint32_t n);