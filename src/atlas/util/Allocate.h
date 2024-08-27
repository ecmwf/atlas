/*
* (C) Copyright 2013 ECMWF.
*
* This software is licensed under the terms of the Apache Licence Version 2.0
* which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
* In applying this licence, ECMWF does not waive the privileges and immunities
* granted to it by virtue of its status as an intergovernmental organisation nor
* does it submit to any jurisdiction.
*/

#pragma once

#include <cstddef>

namespace atlas {
namespace util {

//------------------------------------------------------------------------------

namespace detail {
void allocate_managed(void** ptr, size_t bytes);
void deallocate_managed(void* ptr, size_t bytes);

void allocate_device(void** ptr, size_t bytes);
void deallocate_device(void* ptr, size_t bytes);

void allocate_host(void** ptr, size_t bytes);
void deallocate_host(void* ptr, size_t bytes);

}  // namespace detail

template <typename T>
void allocate_managedmem(T*& data, size_t N) {
    if (N != 0) {
        detail::allocate_managed(reinterpret_cast<void**>(&data), N * sizeof(T));
    }
}

template <typename T>
void delete_managedmem(T*& data, size_t N) {
    if (data) {
        detail::deallocate_managed(data, N * sizeof(T));
        data = nullptr;
    }
}

template <typename T>
void allocate_devicemem(T*& data, size_t N) {
    if (N != 0) {
        detail::allocate_device(reinterpret_cast<void**>(&data), N * sizeof(T));
    }
}

template <typename T>
void delete_devicemem(T*& data, size_t N) {
    if (data) {
        detail::deallocate_device(data, N * sizeof(T));
        data = nullptr;
    }
}

template <typename T>
void allocate_hostmem(T*& data, size_t N) {
    if (N != 0) {
        detail::allocate_host(reinterpret_cast<void**>(&data), N * sizeof(T));
    }
}

template <typename T>
void delete_hostmem(T*& data, size_t N) {
    if (data) {
        detail::deallocate_host(data, N * sizeof(T));
        data = nullptr;
    }
}


//------------------------------------------------------------------------------

extern "C" {
void atlas__allocate_managedmem_double(double*& a, size_t N);
void atlas__allocate_managedmem_float(float*& a, size_t N);
void atlas__allocate_managedmem_int(int*& a, size_t N);
void atlas__allocate_managedmem_long(long*& a, size_t N);
void atlas__deallocate_managedmem_double(double*& a, size_t N);
void atlas__deallocate_managedmem_float(float*& a, size_t N);
void atlas__deallocate_managedmem_int(int*& a, size_t N);
void atlas__deallocate_managedmem_long(long*& a, size_t N);
}

//------------------------------------------------------------------------------

}  // namespace util
}  // namespace atlas
