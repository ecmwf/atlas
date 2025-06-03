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

#include "pluto/memcpy.h"
#include "pluto/stream.h"

namespace pluto {

// ------------------------------------------------------------------------------------------------------------------------------
// Here follow templated copy functions that work with size instead of bytes

template <class T>
void copy_host_to_device(T* device, const T* host, std::size_t size) {
    memcpy_host_to_device(device, host, sizeof(T) * size);
}

template <class T>
void copy_host_to_device(T* device, const T* host, std::size_t size, stream_view s) {
    memcpy_host_to_device(device, host, sizeof(T) * size, s);
}

template <class T>
void copy_host_to_device(T* device, const T* host) {
    copy_host_to_device(device, host, 1);
}

template <class T>
void copy_host_to_device(T* device, const T* host, stream_view s) {
    copy_host_to_device(device, host, 1, s);
}

template <class T>
void copy_host_to_device(T& device, const T& host) {
    copy_host_to_device(&device, &host);
}

template <class T>
void copy_host_to_device(T& device, const T& host, stream_view s) {
    copy_host_to_device(&device, &host, s);
}

template <class T, class D1, class D2>
void copy_host_to_device(std::unique_ptr<T, D1>& device, const std::unique_ptr<T, D2>& host) {
    copy_host_to_device(device.get(), host.get());
}

template <class T, class D1, class D2>
void copy_host_to_device(std::unique_ptr<T, D1>& device, const std::unique_ptr<T, D2>& host, stream_view s) {
    copy_host_to_device(device.get(), host.get(), s);
}

template <class T, class D>
void copy_host_to_device(std::unique_ptr<T, D>& device, const T* host) {
    copy_host_to_device(device.get(), host);
}

template <class T, class D>
void copy_host_to_device(std::unique_ptr<T, D>& device, const T* host, stream_view s) {
    copy_host_to_device(device.get(), host, s);
}

template <class T, class D>
void copy_host_to_device(std::unique_ptr<T, D>& device, const T& host) {
    copy_host_to_device(device.get(), &host);
}

template <class T, class D>
void copy_host_to_device(std::unique_ptr<T, D>& device, const T& host, stream_view s) {
    copy_host_to_device(device.get(), &host, s);
}

template <class T, class D>
void copy_host_to_device(T* device, const std::unique_ptr<T, D>& host) {
    copy_host_to_device(device, host.get());
}

template <class T, class D>
void copy_host_to_device(T* device, const std::unique_ptr<T, D>& host, stream_view s) {
    copy_host_to_device(device, host.get(), s);
}

template <class T, class D>
void copy_host_to_device(T& device, const std::unique_ptr<T, D>& host) {
    copy_host_to_device(&device, host.get());
}

template <class T, class D>
void copy_host_to_device(T& device, const std::unique_ptr<T, D>& host, stream_view s) {
    copy_host_to_device(&device, host.get(), s);
}


template <class T>
void copy_device_to_host(T* host, const T* device, std::size_t size) {
    memcpy_device_to_host(host, device, sizeof(T) * size);
}

template <class T>
void copy_device_to_host(T* host, const T* device, std::size_t size, stream_view s) {
    memcpy_device_to_host(host, device, sizeof(T) * size, s);
}

template <class T>
void copy_device_to_host(T* host, const T* device) {
    copy_device_to_host(host, device, 1);
}

template <class T>
void copy_device_to_host(T* host, const T* device, stream_view s) {
    copy_device_to_host(host, device, 1, s);
}

template <class T>
void copy_device_to_host(T& host, const T& device) {
    copy_device_to_host(&host, &device);
}

template <class T>
void copy_device_to_host(T& host, const T& device, stream_view s) {
    copy_device_to_host(&host, &device, s);
}

template <class T, class D1, class D2>
void copy_device_to_host(std::unique_ptr<T, D1>& host, const std::unique_ptr<T, D2>& device) {
    copy_device_to_host(host.get(), device.get());
}

template <class T, class D1, class D2>
void copy_device_to_host(std::unique_ptr<T, D1>& host, const std::unique_ptr<T, D2>& device, stream_view s) {
    copy_device_to_host(host.get(), device.get(), s);
}

template <class T, class D>
void copy_device_to_host(std::unique_ptr<T, D>& host, const T* device) {
    copy_device_to_host(host.get(), device);
}

template <class T, class D>
void copy_device_to_host(std::unique_ptr<T, D>& host, const T* device, stream_view s) {
    copy_device_to_host(host.get(), device, s);
}

template <class T, class D>
void copy_device_to_host(std::unique_ptr<T, D>& host, const T& device) {
    copy_device_to_host(host.get(), &device);
}

template <class T, class D>
void copy_device_to_host(std::unique_ptr<T, D>& host, const T& device, stream_view s) {
    copy_device_to_host(host.get(), &device, s);
}

template <class T, class D>
void copy_device_to_host(T* host, const std::unique_ptr<T, D>& device) {
    copy_device_to_host(host, device.get());
}

template <class T, class D>
void copy_device_to_host(T* host, const std::unique_ptr<T, D>& device, stream_view s) {
    copy_device_to_host(host, device.get(), s);
}

template <class T, class D>
void copy_device_to_host(T& host, const std::unique_ptr<T, D>& device) {
    copy_device_to_host(&host, device.get());
}

template <class T, class D>
void copy_device_to_host(T& host, const std::unique_ptr<T, D>& device, stream_view s) {
    copy_device_to_host(&host, device.get(), s);
}

template <class T>
void copy_host_to_device_2D(
    T* device, std::size_t device_pitch,
    const T* host, std::size_t host_pitch,
    std::size_t width, std::size_t height) {

    memcpy_host_to_device_2D(
        device, device_pitch * sizeof(T),
        host, host_pitch * sizeof(T),
        width * sizeof(T), height);
}

template <class T>
void copy_device_to_host_2D(
    T* host, std::size_t host_pitch,
    const T* device, std::size_t device_pitch,
    std::size_t width, std::size_t height) {

    memcpy_device_to_host_2D(
        host, host_pitch * sizeof(T),
        device, device_pitch * sizeof(T),
        width * sizeof(T), height);
}

}  // namespace pluto
