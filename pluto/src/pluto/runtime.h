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

#include <cstddef>

#include "pluto/pluto_config.h"

#if HIC_COMPILER && HIC_BACKEND_CUDA
#include <nv/target>
#endif

namespace pluto {
enum class device_type
{
    host,
    device
};

#if HIC_COMPILER
#if defined(__clang__)
__host__ constexpr device_type get_device_type() noexcept {
    return device_type::host;
}

__device__ constexpr device_type get_device_type() noexcept {
    return device_type::device;
}
#else
__host__ __device__ constexpr device_type get_device_type() noexcept {
    NV_IF_TARGET(NV_IS_HOST, (return device_type::host;), (return device_type::device;));
}
#endif

__host__ __device__ constexpr bool is_on_device() noexcept {
    return get_device_type() == device_type::device;
}
__host__ __device__ constexpr bool is_on_host() noexcept {
    return get_device_type() == device_type::host;
}

#else

constexpr device_type get_device_type() noexcept {
    return device_type::host;
}

constexpr bool is_on_device() noexcept {
    return false;
}
constexpr bool is_on_host() noexcept {
    return true;
}

#endif

std::size_t devices();

}  // namespace pluto
