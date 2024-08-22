
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
#include <cstddef>

#include "pluto/pluto.h"

[[maybe_unused]] constexpr std::size_t kb = 1024;
[[maybe_unused]] constexpr std::size_t mb = 1024*kb;
[[maybe_unused]] constexpr std::size_t gb = 1024*mb;

int main(int argc, char* argv[]) {
    std::cout << "BEGIN" << std::endl;

    pluto::TraceOptions::instance().enabled = true;
    bool pinning = true;

    std::unique_ptr<pluto::MemoryPoolResource> memory_pool;
    if (pinning) {
      memory_pool =
        std::make_unique<pluto::MemoryPoolResource>(
        std::make_unique<pluto::TraceMemoryResource>("pinned_upstream",
        std::make_unique<pluto::PinnedMemoryResource>()));
    } else {
      memory_pool =
        std::make_unique<pluto::MemoryPoolResource>(
        std::make_unique<pluto::TraceMemoryResource>("upstream"));
    }

    auto test = [&](std::size_t size, std::size_t n = 1) {
      pluto::allocator<std::byte> allocator{memory_pool.get()};
      std::vector<std::byte*> arrays(n);
      for (std::size_t k=0; k<n; ++k) {
        std::cout << "+ allocate array["<<k<<"], size = " << size/mb << " mb" << std::endl;
        arrays[k] = allocator.allocate(size);
      }
      for (std::size_t k=0; k<n; ++k) {
        std::cout << "- deallocate array["<<n-1-k<<"], size = " << size/mb << " mb" << std::endl;
        allocator.deallocate(arrays[n-1-k],size);
      }
    };

    for (int iter=0; iter<2; ++iter) {
      std::cout << "\n------ iteration " << iter+1 << std::endl;
      test(64*mb, 4);
      test(32*mb);
      test(1*mb);
      test(128*mb);
    }

    std::cout << "END" << std::endl;
}
