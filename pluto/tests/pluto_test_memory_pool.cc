
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
            std::cout << "  allocate " << name_ << " : " << size_ << " bytes" << std::endl;
            data_ = alloc_.allocate(size_);
        }
    }
    ~vector() {
        if (size_) {
            std::cout << "  deallocate " << name_ << " : " << size_ << " bytes" << std::endl;
            alloc_.deallocate(data_, size_);
        }
    }
    std::byte* data_  = nullptr;
    std::size_t size_ = 0;
    pluto::allocator<std::byte> alloc_;
    std::string name_;
};


int main(int argc, char* argv[]) {
    std::cout << "BEGIN" << std::endl;

    pluto::set_trace(true);

    // Uncomment to immediately reserve large chunk of memory
    //   pluto::pool_resource()->reserve(2*gb);

    pluto::allocator<std::byte> allocator(pluto::pool_resource());

    for (int iter = 0; iter < 2; ++iter) {
        vector array1("array1", 200 * mb, allocator);
        vector array2("array2", 200 * mb, allocator);
        vector array3("array3", 1.2 * gb, allocator);
    }


    std::cout << "END" << std::endl;
}
