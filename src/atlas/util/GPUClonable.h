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

#if ATLAS_HAVE_CUDA
#include "hic/hic.h"
#endif

namespace atlas {
namespace util {

template <typename Base>
struct GPUClonable {
    GPUClonable(Base* base_ptr): base_ptr_(base_ptr), gpu_object_ptr_(nullptr) {
#if ATLAS_HAVE_CUDA
        hicMalloc(&gpu_object_ptr_, sizeof(Base));
#endif
    }

    ~GPUClonable() {
        if (gpu_object_ptr_) {
#if ATLAS_HAVE_CUDA
            hicFree(gpu_object_ptr_);
#endif
        }
    }

    Base* gpu_object_ptr() { return static_cast<Base*>(gpu_object_ptr_); }

    void updateDevice() {
#if ATLAS_HAVE_CUDA
        hicMemcpy(gpu_object_ptr_, base_ptr_, sizeof(Base), hicMemcpyHostToDevice);
#endif
    }
    void updateHost() {
#if ATLAS_HAVE_CUDA
        hicMemcpy(base_ptr_, gpu_object_ptr_, sizeof(Base), hicMemcpyDeviceToHost);
#endif
    }

    Base* base_ptr_;
    void* gpu_object_ptr_;
};

}  // namespace util
}  // namespace atlas
