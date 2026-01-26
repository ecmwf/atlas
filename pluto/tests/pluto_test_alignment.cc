
/*
 * (C) Copyright 2024- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "pluto_testing.h"

#include <cstddef>
#include <iostream>
#include <string>
#include <limits>
#include <vector>

#include "pluto/pluto.h"

struct alignas(32) Vec3D_alignment32 {
    double x;
    double y;
    double z;


    static std::string type() { return "Vec3D_alignment32"; }
};

struct Vec3D {
    double x;
    double y;
    double z;

    static std::string type() { return "Vec3D"; }
};

std::size_t get_alignment(void* ptr) {
    for (std::size_t alignment : {2056,1024,512,256,128,64,32,16,8,4}) {
        if (pluto::is_aligned(ptr, alignment)) {
            return alignment;
        }
    }
    return 0;
}


template<typename T>
struct Allocated {
    Allocated(const pluto::allocator<T>& alloc): alloc_(alloc) {
        data = alloc_.allocate(1);
    }
    ~Allocated() {
        alloc_.deallocate(data, 1);
    }
    T* data;
    pluto::allocator<T> alloc_;
};

template<typename Vec>
void run_case() {
    std::cerr << "Case: vector<" << Vec::type() << ">" << std::endl;
    int N = 5;

    std::vector<Vec, pluto::host::allocator<Vec>> vec(N);

    // We expect the pluto::host::allocator to use pluto::host_resource, which aligns to pluto::default_alignment
    EXPECT( pluto::is_aligned(vec.data(), pluto::default_alignment()) );

    std::size_t alignment = std::numeric_limits<std::size_t>::max();
    for (auto& element: vec) {
        EXPECT( pluto::is_aligned(&element, alignof(Vec)) );
        alignment = std::min(alignment, get_alignment(&element));
        std::cerr << "  element alignment = " << get_alignment(&element) << std::endl;
    }
    std::cerr << "  minimum alignment = " << alignment << std::endl;
}

int main([[maybe_unused]] int argc, [[maybe_unused]] char* argv[]) {
    std::cerr << "TEST begin" << std::endl;

    static_assert(alignof(Vec3D_alignment32) == 32);
    run_case<Vec3D>();
    run_case<Vec3D_alignment32>();

    {
        auto test = [](const std::string& alloc_str, const pluto::allocator<int>& allocator, std::size_t expected_alignment) {
            std::cerr << "Case: alignment of " << alloc_str << std::endl;
            Allocated<int> allocated(allocator);
            EXPECT(pluto::is_aligned(allocated.data, expected_alignment));
        };
        #define TEST(allocator, expected_alignment) test(#allocator, allocator, expected_alignment);
        TEST(pluto::allocator<int>(),         alignof(int));               // by default uses pluto::new_delete_resource()
        TEST(pluto::host::allocator<int>(),   pluto::default_alignment()); // by default uses pluto::host_resource()
        TEST(pluto::device::allocator<int>(), pluto::default_alignment()); // by default uses pluto::device_resource()

        TEST(pluto::new_delete_resource(),      alignof(double));
        TEST(pluto::host_resource(),            pluto::default_alignment());
        TEST(pluto::device_resource(),          pluto::default_alignment());
        TEST(pluto::managed_resource(),         pluto::default_alignment());
        TEST(pluto::host_pool_resource(),       pluto::default_alignment());
        TEST(pluto::device_pool_resource(),     pluto::default_alignment());
        TEST(pluto::managed_pool_resource(),    pluto::default_alignment());
    }

    return pluto_testing_return();
}
