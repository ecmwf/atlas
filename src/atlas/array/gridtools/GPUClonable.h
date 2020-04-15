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

namespace atlas {
namespace array {
namespace gridtools {

template <typename Base>
struct GPUClonable {
    GPUClonable( Base* base_ptr ) : base_ptr_( base_ptr ), gpu_object_ptr_( nullptr ) {
#if ATLAS_GRIDTOOLS_STORAGE_BACKEND_CUDA
        cudaMalloc( &gpu_object_ptr_, sizeof( Base ) );
#endif
    }

    ~GPUClonable() {
#if ATLAS_GRIDTOOLS_STORAGE_BACKEND_CUDA
        assert( gpu_object_ptr_ );
        cudaFree( gpu_object_ptr_ );
#endif
    }

    Base* gpu_object_ptr() { return static_cast<Base*>( gpu_object_ptr_ ); }

    void updateDevice() {
#if ATLAS_GRIDTOOLS_STORAGE_BACKEND_CUDA
        cudaMemcpy( gpu_object_ptr_, base_ptr_, sizeof( Base ), cudaMemcpyHostToDevice );
#endif
    }
    void updateHost() {
#if ATLAS_GRIDTOOLS_STORAGE_BACKEND_CUDA
        cudaMemcpy( base_ptr_, gpu_object_ptr_, sizeof( Base ), cudaMemcpyDeviceToHost );
#endif
    }

    Base* base_ptr_;
    void* gpu_object_ptr_;
};

}  // namespace gridtools
}  // namespace array
}  // namespace atlas
