/*
 * (C) Copyright 2024- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <iostream>
#include <algorithm>
#include <vector>
#include <chrono>

#include "pluto/pluto.h"

// ---------------------------------------------------------------------------------------------------
// Helper array types to allocate and deallocate uninitialized memory using allocator.
// This could be in a separate file

template <typename T, typename Allocator>
class array{
public:
  using value_type = T;
  array(std::size_t size) :
    size_{size} {
    data_ = alloc_.allocate(size_);
  }
  template <typename Alloc>
  array(std::size_t size, const Alloc& alloc) :
    alloc_{alloc},
    size_{size} {
    data_ = alloc_.allocate(size_);
  }
  ~array() {
    alloc_.deallocate(data_, size_);
  }
  value_type&       operator[](std::size_t i) { return data_[i]; }
  const value_type& operator[](std::size_t i) const { return data_[i]; }
  value_type*       data() { return data_; };
  const value_type* data() const { return data_; }
  std::size_t size() const { return size_; }
private:
  Allocator alloc_;
  std::size_t size_;
  value_type* data_;
};

template <typename T> using device_array = array<T,pluto::device::allocator<T>>;
template <typename T> using host_array   = array<T,pluto::host::allocator<T>>;

// ---------------------------------------------------------------------------------------------------
// Kernel to add +1 to device array. This could be in a separate CUDA/HIP source file

#if PLUTO_HAVE_HIC
template<typename T>
HIC_GLOBAL
void kernel_plus_one_on_device(T* d, int n) {
    const int idx{int(blockDim.x) * int(blockIdx.x) + int(threadIdx.x)};
    const int stride{int(blockDim.x) * int(gridDim.x)};
    for (int i{idx}; i < n; i += stride) {
      d[i] += 1.;
    }
}

template<typename T>
void plus_one_on_device(T* d, int n, const pluto::Stream& stream) {
    const int threads_per_block{1024};
    const int blocks_per_grid{32};
    kernel_plus_one_on_device<<<blocks_per_grid, threads_per_block, 0, stream.value<hicStream_t>()>>>(d, n);
}
#else
template<typename T>
void plus_one_on_device(T* d, int n, const pluto::Stream& stream) {
  for(int i=0; i<n; ++i) { d[i] += 1.; }
}
#endif


// ---------------------------------------------------------------------------------------------------

int main(int argc, char* argv[]) {
    // Choose data type
    using value_type = float;

    // Total problem size
    std::size_t size = 1'000'000'000;
    if(argc > 1) { size = std::atoi(argv[1]); }
    std::cerr << "size    = " << size << "   (" << size*sizeof(value_type)/1024./1024./1024. << " Gb)" << std::endl;

    // Number of streams to divide problem size in
    std::size_t nb_streams = 7;
    if(argc > 2) { nb_streams = std::atoi(argv[2]); }
    std::cerr << "streams = " << nb_streams << std::endl;

    // A stream pool
    std::vector<pluto::Stream> streams(nb_streams);

    // Set host memory default to be pinned for faster memory transfer to device
    pluto::PinnedMemoryResource pinned_resource;
    pluto::host::set_default_resource(&pinned_resource);

    std::cerr << "host alloc" << std::endl;
    host_array<value_type> array_h1(size);
    host_array<value_type> array_h2(size);

    std::cerr << "device alloc" << std::endl;
    device_array<value_type> array_d1(size);

    std::cerr << "async loop start" << std::endl;
    auto start = std::chrono::steady_clock::now();
    for(std::size_t jstream=0; jstream<streams.size(); ++jstream) {
      pluto::scope stream_scope;
      pluto::set_default_stream(streams[jstream]);

      const auto& stream = streams[jstream];
      const auto& stream_offset = jstream * size/streams.size();
      const auto& stream_size   = (jstream < streams.size()-1 ? size/streams.size() : size - stream_offset);
      auto* h1 = array_h1.data()+stream_offset;
      auto* h2 = array_h2.data()+stream_offset;
      auto* d1 = array_d1.data()+stream_offset;

      device_array<value_type> stream_tmp(stream_size);
      auto* dtmp = stream_tmp.data();

      h1[stream_size-1] = 1.;
      h2[stream_size-1] = -1.;
      pluto::copy_host_to_device(d1, h1, stream_size, stream);
      plus_one_on_device(d1, stream_size, stream);
      pluto::copy_device_to_host(h2, d1, stream_size, stream);
    }
    std::cerr << "async loop end" << std::endl;
    pluto::wait();
    auto end = std::chrono::steady_clock::now();
    std::cout << "execution without allocations took " << std::chrono::duration<double>(end-start).count() << " s" << std::endl;

    if( array_h2[size-1] != 2. ) {
      std::cerr << "ERROR: asynchronous execution not successful, h2[size-1] = " << array_h2[size-1] << std::endl;
    }
}
