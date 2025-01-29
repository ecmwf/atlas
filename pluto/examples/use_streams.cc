/*
 * (C) Copyright 2024- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <algorithm>
#include <chrono>
#include <iostream>
#include <vector>

#include "pluto/pluto.h"

// ---------------------------------------------------------------------------------------------------
// Helper array types to allocate and deallocate uninitialized memory using allocator.
// This could be in a separate file

// SFINAE test
template <typename T>
class has_allocate_async {
    typedef char one;
    struct two {
        char x[2];
    };

    template <typename C>
    static one test(decltype(&C::allocate_async));
    template <typename C>
    static two test(...);

public:
    enum
    {
        value = sizeof(test<T>(0)) == sizeof(char)
    };
};

template <typename T>
class has_deallocate_async {
    typedef char one;
    struct two {
        char x[2];
    };

    template <typename C>
    static one test(decltype(&C::deallocate_async));
    template <typename C>
    static two test(...);

public:
    enum
    {
        value = sizeof(test<T>(0)) == sizeof(char)
    };
};

template <typename T>
static constexpr bool has_async() {
    return has_allocate_async<T>::value && has_deallocate_async<T>::value;
}

template <typename T, typename Allocator = pluto::allocator<T>>
class array {
public:
    using value_type     = T;
    using allocator_type = Allocator;

    array(std::size_t size): size_{size} { data_ = alloc_.allocate(size_); }

    array(std::size_t size, const pluto::stream& stream): size_{size}, stream_(&stream) {
        data_ = alloc_.allocate_async(size_, *stream_);
    }

    template <typename Alloc>
    array(std::size_t size, const Alloc& alloc): alloc_{alloc}, size_{size} {
        data_ = alloc_.allocate(size_);
    }

    template <typename Alloc>
    array(std::size_t size, const pluto::stream& stream, const Alloc& alloc):
        alloc_{alloc}, size_{size}, stream_(&stream) {
        if constexpr (has_async<allocator_type>()) {
            data_ = alloc_.allocate_async(size_, *stream_);
        }
        else {
            data_ = alloc_.allocate(size_);
        }
    }

    ~array() {
        if constexpr (has_async<allocator_type>()) {
            if (stream_) {
                alloc_.deallocate_async(data_, size_, *stream_);
            }
            else {
                alloc_.deallocate(data_, size_);
            }
        }
        else {
            alloc_.deallocate(data_, size_);
        }
    }
    void set_stream(const pluto::stream& stream) { stream_ = &stream; }
    value_type& operator[](std::size_t i) { return data_[i]; }
    const value_type& operator[](std::size_t i) const { return data_[i]; }
    value_type* data() { return data_; };
    const value_type* data() const { return data_; }
    std::size_t size() const { return size_; }

private:
    allocator_type alloc_;
    std::size_t size_;
    value_type* data_{nullptr};
    pluto::stream const* stream_{nullptr};
};

// ---------------------------------------------------------------------------------------------------
// Kernel to add +1 to device array. This could be in a separate CUDA/HIP source file

#if PLUTO_HAVE_HIC
template <typename T>
HIC_GLOBAL void kernel_plus_one_on_device(T* d, int n) {
    const int idx{int(blockDim.x) * int(blockIdx.x) + int(threadIdx.x)};
    const int stride{int(blockDim.x) * int(gridDim.x)};
    for (int i{idx}; i < n; i += stride) {
        d[i] += 1.;
    }
}

template <typename T>
void plus_one_on_device(T* d, int n, const pluto::stream& stream) {
    const int threads_per_block{1024};
    const int blocks_per_grid{32};
    kernel_plus_one_on_device<<<blocks_per_grid, threads_per_block, 0, stream.value<hicStream_t>()>>>(d, n);
}
#else
template <typename T>
void plus_one_on_device(T* d, int n, const pluto::stream& stream) {
    for (int i = 0; i < n; ++i) {
        d[i] += 1.;
    }
}
#endif


// ---------------------------------------------------------------------------------------------------

int main(int argc, char* argv[]) {
    // Choose data type
    using value_type = float;

    // Total problem size
    std::size_t size = 1'000'000'000;
    if (argc > 1) {
        size = std::atoi(argv[1]);
    }
    std::cerr << "size    = " << size << "   (" << size * sizeof(value_type) / 1024. / 1024. / 1024. << " Gb)"
              << std::endl;

    // Number of streams to divide problem size in
    std::size_t nb_streams = 7;
    if (argc > 2) {
        nb_streams = std::atoi(argv[2]);
    }
    std::cerr << "streams = " << nb_streams << std::endl;

    // A stream pool
    std::vector<pluto::stream> streams(nb_streams);

    std::cerr << "host alloc" << std::endl;
    array<value_type> array_h1(size, pluto::pinned_resource());
    array<value_type> array_h2(size, pluto::pinned_resource());

    std::cerr << "device alloc" << std::endl;
    auto& device_resource = *pluto::device_resource();
    array<value_type> array_d1(size, &device_resource);

    std::cerr << "async loop start" << std::endl;
    auto start = std::chrono::steady_clock::now();
    for (std::size_t jstream = 0; jstream < streams.size(); ++jstream) {
        const auto& stream = streams[jstream];
        {
            const auto& stream_offset = jstream * size / streams.size();
            const auto& stream_size   = (jstream < streams.size() - 1 ? size / streams.size() : size - stream_offset);
            array<value_type> stream_tmp(stream_size, stream, &device_resource);

            auto* h1 = array_h1.data() + stream_offset;
            auto* h2 = array_h2.data() + stream_offset;
            auto* d1 = array_d1.data() + stream_offset;

            auto* dtmp = stream_tmp.data();

            h1[stream_size - 1] = 1.;
            h2[stream_size - 1] = -1.;
            pluto::copy_host_to_device(d1, h1, stream_size, stream);
            plus_one_on_device(d1, stream_size, stream);
            pluto::copy_device_to_host(h2, d1, stream_size, stream);
        }
        //stream.wait();
    }
    std::cerr << "async loop end" << std::endl;
    pluto::wait();
    auto end = std::chrono::steady_clock::now();
    std::cout << "execution without allocations took " << std::chrono::duration<double>(end - start).count() << " s"
              << std::endl;

    if (array_h2[size - 1] != 2.) {
        std::cerr << "ERROR: asynchronous execution not successful, h2[size-1] = " << array_h2[size - 1] << std::endl;
    }
}
