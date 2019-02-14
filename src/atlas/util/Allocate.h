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

#include "atlas/library/config.h"

namespace atlas {
namespace util {

//------------------------------------------------------------------------------

namespace detail {
void allocate_cudamanaged( void** ptr, size_t size );
void deallocate_cudamanaged( void* ptr );
}  // namespace detail

template <typename T>
void allocate_managedmem( T*& data, idx_t N ) {
    if ( N != 0 ) {
        detail::allocate_cudamanaged( reinterpret_cast<void**>( &data ), static_cast<size_t>( N ) * sizeof( T ) );
    }
}

template <typename T>
void delete_managedmem( T*& data ) {
    if ( data ) {
        detail::deallocate_cudamanaged( data );
        data = nullptr;
    }
}

//------------------------------------------------------------------------------

extern "C" {
void atlas__allocate_managedmem_double( double*& a, idx_t N );
void atlas__allocate_managedmem_float( float*& a, idx_t N );
void atlas__allocate_managedmem_int( int*& a, idx_t N );
void atlas__allocate_managedmem_long( long*& a, idx_t N );
void atlas__deallocate_managedmem( void*& a );
}

//------------------------------------------------------------------------------

}  // namespace util
}  // namespace atlas
