/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#pragma once

#include "atlas/library/config.h"

#include "pluto/pluto.h"

namespace atlas {
namespace util {

template <typename Base>
struct GPUClonable {
    GPUClonable(Base* base_ptr): base_ptr_(base_ptr), gpu_object_ptr_(nullptr) {
        if constexpr (ATLAS_HAVE_GPU) {
            gpu_object_ptr_ = pluto::device_resource()->allocate(sizeof(Base));
        }
    }

    ~GPUClonable() {
        if (gpu_object_ptr_ != nullptr) {
            pluto::device_resource()->deallocate(gpu_object_ptr_, sizeof(Base));
        }
    }

    Base* gpu_object_ptr() { return static_cast<Base*>(gpu_object_ptr_); }

    void updateDevice() {
        if constexpr (ATLAS_HAVE_GPU) {
            pluto::memcpy_host_to_device(gpu_object_ptr_, base_ptr_, sizeof(Base));
        }
    }
    void updateHost() {
        if constexpr (ATLAS_HAVE_GPU) {
            pluto::memcpy_device_to_host(base_ptr_, gpu_object_ptr_, sizeof(Base));
        }
    }

    Base* base_ptr_;
    void* gpu_object_ptr_;
};

}  // namespace util
}  // namespace atlas
