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

#if ATLAS_STORAGE_BACKEND_CUDA

#include "atlas/library/config.h"

namespace atlas {
namespace array {

template <typename Base>
struct GPUClonable {
    GPUClonable(Base* base_ptr): base_ptr_(base_ptr), gpu_object_ptr_(nullptr) {
        cudaMalloc(&gpu_object_ptr_, sizeof(Base));
    }

    ~GPUClonable() {
        if (gpu_object_ptr_) {
            cudaFree(gpu_object_ptr_);
        }
    }

    Base* gpu_object_ptr() { return static_cast<Base*>(gpu_object_ptr_); }

    void updateDevice() {
        cudaMemcpy(gpu_object_ptr_, base_ptr_, sizeof(Base), cudaMemcpyHostToDevice);
    }
    void updateHost() {
        cudaMemcpy(base_ptr_, gpu_object_ptr_, sizeof(Base), cudaMemcpyDeviceToHost);
    }

    Base* base_ptr_;
    void* gpu_object_ptr_;
};

}  // namespace array
}  // namespace atlas

#endif
