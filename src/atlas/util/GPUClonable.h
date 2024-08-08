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

#include "hic/hic.h"

namespace atlas {
namespace util {

template <typename Base>
struct GPUClonable {
    GPUClonable(Base* base_ptr): base_ptr_(base_ptr), gpu_object_ptr_(nullptr) {
        if constexpr (ATLAS_HAVE_GPU) {
            HIC_CALL(hicMalloc(&gpu_object_ptr_, sizeof(Base)));
        }
    }

    ~GPUClonable() {
        if (gpu_object_ptr_ != nullptr) {
            HIC_CALL(hicFree(gpu_object_ptr_));
        }
    }

    Base* gpu_object_ptr() { return static_cast<Base*>(gpu_object_ptr_); }

    void updateDevice() {
        if constexpr (ATLAS_HAVE_GPU) {
            HIC_CALL(hicMemcpy(gpu_object_ptr_, base_ptr_, sizeof(Base), hicMemcpyHostToDevice));
        }
    }
    void updateHost() {
        if constexpr (ATLAS_HAVE_GPU) {
            HIC_CALL(hicMemcpy(base_ptr_, gpu_object_ptr_, sizeof(Base), hicMemcpyDeviceToHost));
        }
    }

    Base* base_ptr_;
    void* gpu_object_ptr_;
};

}  // namespace util
}  // namespace atlas
