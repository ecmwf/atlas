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

namespace pluto {

bool is_device_accessible(const void* ptr);
bool is_host_accessible(const void* ptr);
bool is_managed(const void* ptr);
bool is_pinned(const void* ptr);
bool is_host(const void* ptr);
bool is_device(const void* ptr);

void* get_registered_device_pointer(const void* host_ptr);
template <typename T>
T* get_registered_device_pointer(const T* host_ptr) {
    return static_cast<T*>(get_registered_device_pointer(static_cast<const void*>(host_ptr)));
}

}  // namespace pluto
