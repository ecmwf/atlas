
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

class vector {
public:
vector(std::size_t n, const pluto::allocator<std::byte>& alloc) :
    size_{n},
    alloc_{alloc} {
    if (size_) {
        data_ = alloc_.allocate(size_);
    }
}
~vector() {
    if(size_) {
        alloc_.deallocate(data_,size_);
    }
}
std::byte* data_ = nullptr;
std::size_t size_ = 0;
pluto::allocator<std::byte> alloc_;
};


int main(int argc, char* argv[]) {
    std::cout << "BEGIN" << std::endl;

    pluto::TraceOptions::instance().enabled = true;

    pluto::allocator<std::byte> allocator(pluto::pool_resource());
    // using vector = pluto::host::vector<std::byte>;

    vector array1(200*mb,allocator);
    vector array2(200*mb,allocator);
    vector array3(1.2*gb,allocator);


    std::cout << "END" << std::endl;
}
