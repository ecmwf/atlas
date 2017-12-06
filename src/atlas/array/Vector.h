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
class Vector {
public:
    Vector(size_t N = 0) : data_(N ? new T[N] : nullptr), data_gpu_(nullptr), size_(N) {}

  void resize_impl(size_t N) {
      if( data_gpu_ ) throw eckit::AssertionFailed("we can not resize a vector after has been cloned to device");
      assert(N >= size_);
      if (N == size_) return;

      T* d_ = new T[N];
      for(unsigned int c=0; c < size_; ++c) {
          d_[c] = data_[c];
      }
      if( data_ ) delete[] data_;
      data_ = d_;
  }

  void resize(size_t N) {
    resize_impl(N);
    size_ = N;
  }

  void resize(size_t N, T val) {
    resize_impl(N);
    for(unsigned int c=size_; c < N; ++c) {
        data_[c] = val;
    }

    size_ = N;
  }

  void cloneToDevice() {
    if(!data_gpu_) {
#if ATLAS_GRIDTOOLS_STORAGE_BACKEND_CUDA
        ::cudaMalloc((void**)(&data_gpu_), sizeof(T*)*size_);

        T* buff = new T[size_];

        for(size_t i=0; i < size(); ++i) {
          data_[i]->cloneToDevice();
          buff[i] = data_[i]->gpu_object_ptr();
        }
        ::cudaMemcpy(data_gpu_, buff, sizeof(T*)*size_, cudaMemcpyHostToDevice);
        delete buff;
#else
        data_gpu_ = data_;
#endif
        size_gpu_ = size_;
    }
    else {
        assert(size_gpu_ == size_);
#if ATLAS_GRIDTOOLS_STORAGE_BACKEND_CUDA
        for(size_t i=0; i < size(); ++i) {
          data_[i]->cloneToDevice();
          assert(data_gpu_[i] == data_[i]->gpu_object_ptr());
        }
#endif
    }
  }
  void cloneFromDevice() {
    assert(data_gpu_);

#if ATLAS_GRIDTOOLS_STORAGE_BACKEND_CUDA

    for(size_t i=0; i < size(); ++i) {
      data_[i]->cloneFromDevice();
    }
#endif
  }

  T* gpu_object_ptr() { return data_gpu_;}

  T* data() { return data_;}
  T* data_gpu() { return data_gpu_;}

  size_t size() const { return size_;}

private:
  T* data_;
  T* data_gpu_;
  size_t size_;
  size_t size_gpu_;
};

template <typename T>
class VectorView {
public:
  VectorView() : vector_(NULL),data_(NULL), size_(0) {}
  VectorView(Vector<T> const& vector, T* data) : vector_(&vector), data_(data), size_(vector.size()) {}

  ATLAS_HOST_DEVICE
  T& operator[](size_t idx) {
      assert(idx < size_);

      return data_[idx];
  }

  ATLAS_HOST_DEVICE
  T const& operator[](size_t idx) const {
      assert(idx < size_);
      return data_[idx];
  }
  ATLAS_HOST_DEVICE
  T base() { return *data_; }

  ATLAS_HOST_DEVICE
  size_t size() const { return size_;}

  bool is_valid(Vector<T>& vector) {
    return (&vector) == vector_ && (data_ != NULL);
  }
public:
  Vector<T> const* vector_;
  T* data_;
  size_t size_;
};

template <typename T>
VectorView<T>
make_host_vector_view(Vector<T> vector_) {
  return VectorView<T>(vector_, vector_.data());
}

template <typename T>
VectorView<T>
make_device_vector_view(Vector<T> vector_) {
  return VectorView<T>(vector_, vector_.data_gpu());
}

//------------------------------------------------------------------------------

} // namespace array
} // namespace atlas
