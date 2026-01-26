
/*
 * (C) Copyright 2024- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <cstddef>
#include <iostream>

#include "pluto/pluto.h"

[[maybe_unused]] constexpr std::size_t kb = 1024;
[[maybe_unused]] constexpr std::size_t mb = 1024 * kb;
[[maybe_unused]] constexpr std::size_t gb = 1024 * mb;

int main([[maybe_unused]] int argc, [[maybe_unused]] char* argv[]) {
    std::cout << "BEGIN" << std::endl;

    pluto::trace::enable(true);
    // This could also have been done with environment variable:
    //     export PLUTO_TRACE=1

    // The default pool options can be changed before using the memory pool for the first time.
    // If not changed, the default configuration allows for *all* allocations to be handled by a memory pool,
    // no matter the size.
    pluto::set_default_pool_options({
        /* .max_blocks_per_chunk = */ 4,
        /* .largest_required_pool_block = */ 64*mb
    });
    // This modification lets only allocations <= 64M be handled by the memory pool.
    // When the pool needs to grow, 4 chunks of pool memory are allocated at a time.
    // This could also have been achieved with environment variables:
    //     export PLUTO_MAX_BLOCKS_PER_CHUNK=4
    //     export PLUTO_LARGEST_REQUIRED_POOL_BLOCK=64M

    pluto::host::set_default_resource( pluto::pinned_pool_resource() );
    // This can also be set with environment variable:
    //     export PLUTO_HOST_DEFAULT_MEMORY_RESOURCE=pluto::host_pool_resource

    auto test = [&](std::size_t size, std::size_t n = 1) {
        auto var_name = [](int k) {
            return std::string("var_") + std::to_string(k);
        };
        pluto::host::allocator<std::byte> allocator;
        std::vector<std::byte*> var_data(n);
        for (std::size_t k = 0; k < n; ++k) {
            var_data[k] = allocator.allocate(var_name(k),size);
        }
        for (std::size_t k = n; k-- > 0; ) {
            allocator.deallocate(var_name(k), var_data[k], size);
        }
    };

    for (int iter = 0; iter < 2; ++iter) {
        std::cout << "\n------ iteration " << iter + 1 << std::endl;
        test(64 * mb, 4);
        test(32 * mb);
        test(1 * mb);
        test(128 * mb);
    }

    std::cout << "\nEND\n" << std::endl;

    // Release all memory pools
    pluto::release();

    std::cout << '\n' << pluto::memory::report() << std::endl;

}
