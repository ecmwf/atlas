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

#include <memory>

#include "pluto/offload/copy.h"
#include "pluto/host/allocator.h"
#include "pluto/host/unique_ptr.h"

namespace pluto::host {

template <class T> 
auto make_copy(const T* device) {
    allocator<T> alloc;
    T* p = alloc.allocate(1);
    copy_device_to_host(p, device);
    return host::unique_ptr<T>(p, alloc);
}
template <class T, class D> 
auto make_copy(const std::unique_ptr<T,D>& device) {
    return make_copy(device.get());
}
template <class T> 
auto make_copy(const T& device) {
    return make_copy(&device);
}

}
