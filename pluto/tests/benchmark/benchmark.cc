#include "benchmark.h"

#include <chrono>
#include <cstring>
#include <iostream>

float launch_benchmark(float* d_output, float const* d_input_1, float const* d_input_2, uint32_t n) {
    const auto start{std::chrono::steady_clock::now()};
    launch_benchmark_kernel(n, [=] HIC_DEVICE(uint32_t i) { d_output[i] = d_input_1[i] + d_input_2[i]; });
    const auto end{std::chrono::steady_clock::now()};
    return std::chrono::duration<float>{end - start}.count();
}

float launch_benchmark(double* d_output, double const* d_input_1, double const* d_input_2, uint32_t n) {
    const auto start{std::chrono::steady_clock::now()};
    launch_benchmark_kernel(n, [=] HIC_DEVICE(uint32_t i) { d_output[i] = d_input_1[i] + d_input_2[i]; });
    const auto end{std::chrono::steady_clock::now()};
    return std::chrono::duration<float>{end - start}.count();
}


float touch_host_memory(float* h_buffer, uint32_t n) {
    const auto start{std::chrono::steady_clock::now()};
    std::memset(h_buffer, 0, n * sizeof(float));
    const auto end{std::chrono::steady_clock::now()};
    return std::chrono::duration<float>{end - start}.count();
}

float touch_host_memory(double* h_buffer, uint32_t n) {
    const auto start{std::chrono::steady_clock::now()};
    std::memset(h_buffer, 0, n * sizeof(double));
    const auto end{std::chrono::steady_clock::now()};
    return std::chrono::duration<float>{end - start}.count();
}


float run_timed(std::function<void()> function, int num_repeats, int num_warmups) {
    for (int i{0}; i < num_warmups; ++i) {
        function();
    }

    HIC_CALL(hicDeviceSynchronize());

    const auto cstart{std::chrono::steady_clock::now()};
    for (int i{0}; i < num_repeats; ++i) {
        function();
    }
    const auto cend{std::chrono::steady_clock::now()};

    HIC_CALL(hicDeviceSynchronize());
    return std::chrono::duration<float>{cend - cstart}.count() / float(num_repeats);
}
