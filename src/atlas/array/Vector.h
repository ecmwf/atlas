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

namespace atlas {
namespace array {

//------------------------------------------------------------------------------

template <typename T>
class Vector {
public:
  Vector() : data_(NULL), size_(0) {}
  Vector(size_t N) : data_(new T[N]()), size_(N) {}

  void resize(size_t N) {
    assert(N >= size_);
    if (N == size_) return;

    T* d_ = new T[N]();
    for(unsigned int c=0; c < size_; ++c) {
        d_[c] = data_[c];
    }
    delete data_;
    data_ = d_;
    size_ = N;
  }

  T& operator[](size_t idx) {
      assert(idx < size_);
      return data_[idx];
  }

  T const& operator[](size_t idx) const {
      assert(idx < size_);
      return data_[idx];
  }

  size_t size() { return size_;}

private:
  T* data_;
  size_t size_;
};

//------------------------------------------------------------------------------

} // namespace array
} // namespace atlas
