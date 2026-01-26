
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
#include <string>
#include <string_view>

#include "pluto/pluto.h"

[[maybe_unused]] constexpr std::size_t kb = 1024;
[[maybe_unused]] constexpr std::size_t mb = 1024 * kb;
[[maybe_unused]] constexpr std::size_t gb = 1024 * mb;

class vector {
public:
    vector(std::string_view name, std::size_t n, const pluto::allocator<std::byte>& alloc):
        name_(name), size_{n}, alloc_{alloc} {
        if (size_) {
            data_ = alloc_.allocate(name_, size_);
        }
    }
    ~vector() {
        if (size_) {
            alloc_.deallocate(name_, data_, size_);
        }
    }
    std::byte* data_ = nullptr;
    std::string name_;
    std::size_t size_ = 0;
    pluto::allocator<std::byte> alloc_;
};


int main([[maybe_unused]] int argc, [[maybe_unused]] char* argv[]) {
    std::cout << "BEGIN" << std::endl;

    pluto::trace::enable(true);

    // Uncomment to immediately reserve large chunk of memory
    //   pluto::host_pool_resource()->reserve(2*gb);

    auto memory_pool = pluto::TraceMemoryResource(pluto::host_pool_resource());
    pluto::set_default_resource(&memory_pool);

    pluto::allocator<std::byte> allocator;

    for (int iter = 0; iter < 2; ++iter) {
        vector array1("array1", 200 * mb, allocator);
        vector array2("array2", 200 * mb, allocator);
        vector array3("array3", 1.2 * gb, allocator);
    }


    std::cout << "END" << std::endl;
}
