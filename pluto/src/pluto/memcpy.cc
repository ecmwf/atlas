
/*
 * (C) Copyright 2024- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "memcpy.h"

#include <cstring>  // std::memcpy
#include <exception>
#include <iostream>

#include "hic/hic.h"

#include "pluto/pluto_config.h"
#include "pluto/runtime.h"
#include "pluto/wait.h"

#define LOG PLUTO_DEBUGGING

namespace pluto {

static void throw_not_implemented(const char* message, const char* file, int line) {
    std::stringstream err;
    err << "ERROR (NOT IMPLEMENTED) : " << message << " [ at " << file << " : " << line << "]";
    throw std::runtime_error(err.str());
}
#define THROW_NOT_IMPLEMENTED(message) throw_not_implemented(message, __FILE__, __LINE__)

void memcpy_host_to_device(void* device_ptr, const void* host_ptr, std::size_t bytes) {
    if (device_ptr == host_ptr) {
        if (LOG) {
            std::cout << "               > hicMemcpyHostToDevice(device_ptr, host_ptr, bytes:" << bytes
                      << ") --> SKIP because same pointers" << std::endl;
        }
        wait();
        return;
    }
    if (LOG) {
        std::cout << "               > hicMemcpyHostToDevice(device_ptr, host_ptr, bytes:" << bytes << ")" << std::endl;
    }
    if constexpr (PLUTO_HAVE_HIC) {
        HIC_CALL(hicMemcpy(device_ptr, host_ptr, bytes, hicMemcpyHostToDevice));
    }
    else {
        std::memcpy(device_ptr, host_ptr, bytes);
    }
}

void memcpy_host_to_device(void* device_ptr, const void* host_ptr, std::size_t bytes, stream_view s) {
    if (device_ptr == host_ptr) {
        if (LOG) {
            std::cout << "               > hicMemcpyHostToDeviceAsync(device_ptr, host_ptr, bytes:" << bytes
                      << ", stream:" << s.value() << ") --> SKIP because same pointers" << std::endl;
        }
        return;
    }
    if (LOG) {
        std::cout << "               > hicMemcpyHostToDeviceAsync(device_ptr, host_ptr, bytes:" << bytes
                  << ", stream:" << s.value() << ")" << std::endl;
    }
    if constexpr (PLUTO_HAVE_HIC) {
        if (devices()) {
            HIC_CALL(hicMemcpyAsync(device_ptr, host_ptr, bytes, hicMemcpyHostToDevice, s.value<hicStream_t>()));
        }
        else {
            std::memcpy(device_ptr, host_ptr, bytes);
        }
    }
    else {
        std::memcpy(device_ptr, host_ptr, bytes);
    }
}

void memcpy_device_to_host(void* host_ptr, const void* device_ptr, std::size_t bytes) {
    if (device_ptr == host_ptr) {
        if (LOG) {
            std::cout << "               < hicMemcpyDeviceToHost(host_ptr, device_ptr, bytes:" << bytes
                      << ") --> SKIP because same pointers" << std::endl;
        }
        wait();
        return;
    }
    if (LOG) {
        std::cout << "               < hicMemcpyDeviceToHost(host_ptr, device_ptr, bytes:" << bytes << ")" << std::endl;
    }
    if constexpr (PLUTO_HAVE_HIC) {
        HIC_CALL(hicMemcpy(host_ptr, device_ptr, bytes, hicMemcpyDeviceToHost));
    }
    else {
        std::memcpy(host_ptr, device_ptr, bytes);
    }
}

void memcpy_device_to_host(void* host_ptr, const void* device_ptr, std::size_t bytes, stream_view s) {
    if (device_ptr == host_ptr) {
        if (LOG) {
            std::cout << "               < hicMemcpyDeviceToHostAsync(host_ptr, device_ptr, bytes:" << bytes
                      << ", stream:" << s.value<hicStream_t>() << ") --> SKIP because same pointers" << std::endl;
        }
        return;
    }
    if (LOG) {
        std::cout << "               < hicMemcpyDeviceToHostAsync(host_ptr, device_ptr, bytes:" << bytes
                  << ", stream:" << s.value<hicStream_t>() << ")" << std::endl;
    }
    if constexpr (PLUTO_HAVE_HIC) {
        if (devices()) {
            HIC_CALL(hicMemcpyAsync(host_ptr, device_ptr, bytes, hicMemcpyDeviceToHost, s.value<hicStream_t>()));
        }
        else {
            std::memcpy(host_ptr, device_ptr, bytes);
        }
    }
    else {
        std::memcpy(host_ptr, device_ptr, bytes);
    }
}

void memcpy_host_to_device_2D(void* device_ptr,
                              [[maybe_unused]] std::size_t device_pitch_bytes /*stride in bytes to next contiguous chunk on device*/,
                              const void* host_ptr,
                              [[maybe_unused]] std::size_t host_pitch_bytes /*stride in bytes to next contiguous chunk on host*/,
                              [[maybe_unused]] std::size_t width_bytes /*bytes of contiguous chunk*/,
                              [[maybe_unused]] std::size_t height_count /*count of contiguous chunks*/) {
    if (device_ptr == host_ptr) {
        if (LOG) {
            std::cout << "               > hicMemcpyHostToDevice2D(host_ptr, device_ptr) --> SKIP because same pointers"
                      << std::endl;
        }
        wait();
        return;
    }
    if (LOG) {
        std::cout << "               > hicMemcpyHostToDevice2D(host_ptr, device_ptr)" << std::endl;
    }
    if constexpr (PLUTO_HAVE_HIC) {
        HIC_CALL(hicMemcpy2D(device_ptr, device_pitch_bytes, host_ptr, host_pitch_bytes, width_bytes, height_count,
                             hicMemcpyHostToDevice));
    }
    else {
        THROW_NOT_IMPLEMENTED("memcpy_host_to_device_2D for CPU backend");
    }
}

void memcpy_host_to_device_2D(void* device_ptr,
                              [[maybe_unused]] std::size_t device_pitch_bytes /*stride in bytes to next contiguous chunk on device*/,
                              const void* host_ptr,
                              [[maybe_unused]] std::size_t host_pitch_bytes /*stride in bytes to next contiguous chunk on host*/,
                              [[maybe_unused]] std::size_t width_bytes /*bytes of contiguous chunk*/,
                              [[maybe_unused]] std::size_t height_count /*count of contiguous chunks*/,
                              [[maybe_unused]] stream_view s) {
    if (device_ptr == host_ptr) {
        if (LOG) {
            std::cout
                << "               > hicMemcpyHostToDevice2DAsync(host_ptr, device_ptr) --> SKIP because same pointers"
                << std::endl;
        }
        return;
    }
    if (LOG) {
        std::cout << "               > hicMemcpyHostToDevice2DAsync(host_ptr, device_ptr)" << std::endl;
    }
    if constexpr (PLUTO_HAVE_HIC) {
        HIC_CALL(hicMemcpy2DAsync(device_ptr, device_pitch_bytes, host_ptr, host_pitch_bytes, width_bytes, height_count,
                                  hicMemcpyHostToDevice, s.value<hicStream_t>()));
    }
    else {
        THROW_NOT_IMPLEMENTED("memcpy_host_to_device_2D for CPU backend");
    }
}

void memcpy_device_to_host_2D(void* host_ptr,
                              [[maybe_unused]] std::size_t host_pitch_bytes /*stride in bytes to next contiguous chunk on host*/,
                              const void* device_ptr,
                              [[maybe_unused]] std::size_t device_pitch_bytes /*stride in bytes to next contiguous chunk on device*/,
                              [[maybe_unused]] std::size_t width_bytes /*bytes of contiguous chunk*/,
                              [[maybe_unused]] std::size_t height_count /*count of contiguous chunks*/) {
    if (device_ptr == host_ptr) {
        if (LOG) {
            std::cout << "               < hicMemcpyDeviceToHost2D(host_ptr, device_ptr) --> SKIP because same pointers"
                      << std::endl;
        }
        wait();
        return;
    }
    if (LOG) {
        std::cout << "               < hicMemcpyDeviceToHost2D(host_ptr, device_ptr)" << std::endl;
    }

    if constexpr (PLUTO_HAVE_HIC) {
        HIC_CALL(hicMemcpy2D(host_ptr, host_pitch_bytes, device_ptr, device_pitch_bytes, width_bytes, height_count,
                             hicMemcpyDeviceToHost));
    }
    else {
        THROW_NOT_IMPLEMENTED("memcpy_device_to_host_2D for CPU backend");
    }
}
void memcpy_device_to_host_2D(void* host_ptr,
                              [[maybe_unused]] std::size_t host_pitch_bytes /*stride in bytes to next contiguous chunk on host*/,
                              const void* device_ptr,
                              [[maybe_unused]] std::size_t device_pitch_bytes /*stride in bytes to next contiguous chunk on device*/,
                              [[maybe_unused]] std::size_t width_bytes /*bytes of contiguous chunk*/,
                              [[maybe_unused]] std::size_t height_count /*count of contiguous chunks*/,
                              [[maybe_unused]] stream_view s) {
    if (device_ptr == host_ptr) {
        if (LOG) {
            std::cout
                << "               < hicMemcpyDeviceToHost2DAsync(host_ptr, device_ptr) --> SKIP because same pointers"
                << std::endl;
        }
        return;
    }
    if (LOG) {
        std::cout << "               < hicMemcpyDeviceToHost2DAsync(host_ptr, device_ptr)" << std::endl;
    }
    if constexpr (PLUTO_HAVE_HIC) {
        HIC_CALL(hicMemcpy2DAsync(host_ptr, host_pitch_bytes, device_ptr, device_pitch_bytes, width_bytes, height_count,
                                  hicMemcpyDeviceToHost, s.value<hicStream_t>()));
    }
    else {
        THROW_NOT_IMPLEMENTED("memcpy_device_to_host_2D for CPU backend");
    }
}

}  // namespace pluto
