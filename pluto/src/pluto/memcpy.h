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
#include <memory>

#include "pluto/stream.h"

namespace pluto {

void memcpy_host_to_device(void* device_ptr, const void* host_ptr, std::size_t bytes);
void memcpy_host_to_device(void* device_ptr, const void* host_ptr, std::size_t bytes, stream_view);

void memcpy_device_to_host(void* host_ptr, const void* device_ptr, std::size_t bytes);
void memcpy_device_to_host(void* host_ptr, const void* device_ptr, std::size_t bytes, stream_view);

void memcpy_host_to_device_2D(void* device_ptr,
                              std::size_t device_pitch_bytes /*stride in bytes to next contiguous chunk on device*/,
                              const void* host_ptr,
                              std::size_t host_pitch_bytes /*stride in bytes to next contiguous chunk on host*/,
                              std::size_t width_bytes /*bytes of contiguous chunk*/,
                              std::size_t height_count /*count of contiguous chunks*/);
void memcpy_host_to_device_2D(void* device_ptr,
                              std::size_t device_pitch_bytes /*stride in bytes to next contiguous chunk on device*/,
                              const void* host_ptr,
                              std::size_t host_pitch_bytes /*stride in bytes to next contiguous chunk on host*/,
                              std::size_t width_bytes /*bytes of contiguous chunk*/,
                              std::size_t height_count /*count of contiguous chunks*/, stream_view);

void memcpy_device_to_host_2D(void* host_ptr,
                              std::size_t host_pitch_bytes /*stride in bytes to next contiguous chunk on host*/,
                              const void* device_ptr,
                              std::size_t device_pitch_bytes /*stride in bytes to next contiguous chunk on device*/,
                              std::size_t width_bytes /*bytes of contiguous chunk*/,
                              std::size_t height_count /*count of contiguous chunks*/);
void memcpy_device_to_host_2D(void* host_ptr,
                              std::size_t host_pitch_bytes /*stride in bytes to next contiguous chunk on host*/,
                              const void* device_ptr,
                              std::size_t device_pitch_bytes /*stride in bytes to next contiguous chunk on device*/,
                              std::size_t width_bytes /*bytes of contiguous chunk*/,
                              std::size_t height_count /*count of contiguous chunks*/, stream_view);

}  // namespace pluto
