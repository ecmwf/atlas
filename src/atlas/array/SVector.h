/*
* (C) Copyright 1996-2016 ECMWF.
*
* This software is licensed under the terms of the Apache Licence Version 2.0
* which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
* In applying this licence, ECMWF does not waive the privileges and immunities
* granted to it by virtue of its status as an intergovernmental organisation nor
* does it submit to any jurisdiction.
*/

#pragma once

#include <cstddef>
#include <cassert>

#include "atlas/library/config.h"

#if ATLAS_GRIDTOOLS_STORAGE_BACKEND_CUDA
#include <cuda_runtime.h>
#endif

#include "atlas/runtime/ErrorHandling.h"
#include "atlas/library/config.h"

namespace atlas {
namespace array {

//------------------------------------------------------------------------------

template <typename T>
class SVector {
public:
  SVector() : data_(nullptr), size_(0), externally_allocated_(false) {}

  ATLAS_HOST_DEVICE
  SVector(SVector const & other) : data_(other.data_), size_(other.size_), externally_allocated_(other.externally_allocated_){}

  ATLAS_HOST_DEVICE
  SVector(T* data, size_t size): data_(data), size_(size) {}

  SVector(size_t N) : size_(N) {
#if ATLAS_GRIDTOOLS_STORAGE_BACKEND_CUDA
      cudaError_t err = cudaMallocManaged(&data_, N * sizeof(T));
      if(err != cudaSuccess)
          throw eckit::AssertionFailed("failed to allocate GPU memory");
#else
    data_ = (T*)malloc(N*sizeof(T));
#endif
  }
  ATLAS_HOST_DEVICE
  ~SVector(){
#ifndef __CUDA_ARCH__
    if(!externally_allocated_)
      delete_managedmem(data_);
#endif
  }

  void delete_managedmem(T*& data) {
    if( data )  {
#if ATLAS_GRIDTOOLS_STORAGE_BACKEND_CUDA
      cudaError_t err = cudaDeviceSynchronize();
      if(err != cudaSuccess)
          throw eckit::AssertionFailed("failed to synchronize memory");

      err = cudaFree(data);
// The following throws an invalid device memory

      if(err != cudaSuccess)
          throw eckit::AssertionFailed("failed to free GPU memory");

#else
      free(data_);
#endif
      data_ = NULL;
    }
  }
  T* data() { return data_; }
  T const * data() const { return data_; }

  ATLAS_HOST_DEVICE
  T& operator()(const size_t idx) { 
    assert(data_ && idx < size_);
    return data_[idx]; 
  }
  ATLAS_HOST_DEVICE
  T const& operator()(const size_t idx) const { 
    assert(data_ && idx < size_);
    return data_[idx]; 
  }

  ATLAS_HOST_DEVICE
  T& operator[](const size_t idx) { 
    assert(data_ && idx < size_);
    return data_[idx]; 
  }
  ATLAS_HOST_DEVICE
  T const& operator[](const size_t idx) const { 
    assert(data_ && idx < size_);
    return data_[idx]; 
  }

  ATLAS_HOST_DEVICE
  size_t size() const { return size_; }

  void resize_impl(size_t N) {
      assert(N >= size_);
      if (N == size_) return;

      T* d_;
#if ATLAS_GRIDTOOLS_STORAGE_BACKEND_CUDA
      cudaError_t err = cudaMallocManaged(&d_, sizeof(T)*N );
      if(err != cudaSuccess)
          throw eckit::AssertionFailed("failed to allocate GPU memory");

#else
      d_ = (T*)malloc(sizeof(T)*N );
#endif
      for(unsigned int c=0; c < size_; ++c) {
         d_[c] = data_[c];
      }
      delete_managedmem(data_);
      data_ = d_;
  }

  void resize(size_t N) {
    resize_impl(N);
    size_ = N;
  }

private:
  T* data_;
  size_t size_;
  bool externally_allocated_;
};

//------------------------------------------------------------------------------

}  // namespace array
}  // namespace atlas
