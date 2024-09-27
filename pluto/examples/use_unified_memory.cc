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
#include <string_view>

#include "pluto/pluto.h"
#include "hic/hic.h"

// ---------------------------------------------------------------------------------------------------
// Kernel to add +1 to device array. This could be in a separate CUDA/HIP source file

template<typename T>
void plus_one_on_host(T* d, int n) {
  std::cout << "               ~ plus_one_on_host" << std::endl;
  for(int i=0; i<n; ++i) { d[i] += 1.; }
}

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
void plus_one_on_device(T* d, int n) {
    std::cout << "               ~ plus_one_on_device" << std::endl;
    const int threads_per_block{1024};
    const int blocks_per_grid{32};
    HIC_CHECK_KERNEL_LAUNCH();
    kernel_plus_one_on_device<<<blocks_per_grid, threads_per_block>>>(d, n);
    // HIC_CHECK_KERNEL_LAUNCH();
}
#else
template<typename T>
void plus_one_on_device(T* d, int n) {
  plus_one_on_host(d, n);
}
#endif

// ---------------------------------------------------------------------------------------------------

template<typename T>
void print_array(std::string_view symbol, T* d, int n) {
    std::cout << '\n' << symbol << "[0..size] = ";
    for (std::size_t i=0; i<5; ++i) {
      std::cout << d[i] << " ";
    }
    std::cout << "... ";
    for (std::size_t i=n-2; i<n; ++i) {
      std::cout << d[i] << " ";
    }
    std::cout << '\n' << std::endl;
}

void print_info(const void* ptr, bool debug = false) {

    std::cout << "Pointer introspection [" << ptr << "]" << std::endl;
    std::cout << "  pluto::is_managed:           " << pluto::is_managed(ptr) << std::endl;
    std::cout << "  pluto::is_host:              " << pluto::is_host(ptr) << std::endl;
    std::cout << "  pluto::is_device:            " << pluto::is_device(ptr) << std::endl;
    std::cout << "  pluto::is_pinned:            " << pluto::is_pinned(ptr) << std::endl;
    std::cout << "  pluto::is_device_accessible: " << pluto::is_device_accessible(ptr) << std::endl;
    std::cout << "  pluto::is_host_accessible:   " << pluto::is_host_accessible(ptr) << std::endl;

    if (debug) {

#if HIC_BACKEND_HIP && HIP_VERSION_MAJOR < 6
    {
      // Debugging HIP backend pointer attributes because it did not map one-to-one with CUDA
      hipPointerAttribute_t attr_;
      auto err = hipPointerGetAttributes(&attr_, ptr);
      if( err != hipSuccess ) {
        std::cout << "Warning: could not access hipPointerGetAttributes" << std::endl;
        return;
      }

      std::cout << "\nhipPointerAttribute_t: " << std::endl;
      std::cout << "  memoryType:    " << attr_.memoryType << std::endl;
      std::cout << "  device:        " << attr_.device << std::endl;
      std::cout << "  devicePointer: " << attr_.devicePointer << std::endl;
      std::cout << "  hostPointer:   " << attr_.hostPointer << std::endl;
      std::cout << "  isManaged:     " << attr_.isManaged << std::endl;

      std::cout<< "Legend(hipMemoryType):" << std::endl;
      std::cout << "   hipMemoryTypeHost:    " << hipMemoryTypeHost << std::endl;
      std::cout << "   hipMemoryTypeDevice:  " << hipMemoryTypeDevice << std::endl;
      std::cout << "   hipMemoryTypeUnified: " << hipMemoryTypeUnified << std::endl;
    }
#endif

#if PLUTO_HAVE_HIC
    {
      // Debugging HIC
      hicPointerAttributes attr_;
      auto err = hicPointerGetAttributes(&attr_, ptr);
      if( err != hicSuccess ) {
        std::cout << "Warning: could not access hicPointerGetAttributes" << std::endl;
        return;
      }

      std::cout << "\nhicPointerAttributes: " << std::endl;
      std::cout << "  type:          " << attr_.type << std::endl;
      std::cout << "  device:        " << attr_.device << std::endl;
      std::cout << "  devicePointer: " << attr_.devicePointer << std::endl;
      std::cout << "  hostPointer:   " << attr_.hostPointer << std::endl;

      std::cout<< "Legend(hicMemoryType):" << std::endl;
      std::cout << "   hicMemoryTypeUnregistered:  " << hicMemoryTypeUnregistered << std::endl;
      std::cout << "   hicMemoryTypeHost:          " << hicMemoryTypeHost << std::endl;
      std::cout << "   hicMemoryTypeDevice:        " << hicMemoryTypeDevice << std::endl;
      std::cout << "   hicMemoryTypeManaged:       " << hicMemoryTypeManaged << std::endl;
    }
#endif

    }
}

// ---------------------------------------------------------------------------------------------------

enum class Memory {
  Managed,
  Pinned
};

int main(int argc, char* argv[]) {
    // Choose data type
    using value_type = float;

    // Total problem size
    std::size_t size = 1'000;

    // Choose Memory allocation strategy
    Memory strategy = Memory::Managed;

    std::cout << "pluto::devices: " << pluto::devices() << std::endl;

    pluto::host::allocator<value_type> allocator{
      strategy == Memory::Managed ? pluto::managed_resource() :
      strategy == Memory::Pinned  ? pluto::pinned_resource() :
      nullptr
    };

    value_type* h_array = allocator.allocate(size);
    value_type* d_array = pluto::get_registered_device_pointer(h_array);

    print_info(h_array);
    print_info(d_array);

    // Fill array with zero
    if (pluto::is_host_accessible(h_array)) {
      std::fill(h_array,h_array+size,0.);
    }

    print_array("array", h_array, size);

    // Add +1 to each element of array
    if (pluto::is_device_accessible(d_array)) {
      plus_one_on_device(d_array, size);
      pluto::wait();
    }
    else {
      plus_one_on_host(d_array,size);
    }

    print_array("array", h_array, size);

    allocator.deallocate(h_array, size);
    // No need to deallocate d_array !
}
