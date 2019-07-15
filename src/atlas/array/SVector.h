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

#include <algorithm>
#include <cassert>
#include <cstddef>

#include "atlas/library/config.h"
#include "atlas/util/Allocate.h"

namespace atlas {
namespace array {

//------------------------------------------------------------------------------

template <typename T>
class SVector {
public:
    ATLAS_HOST_DEVICE
    SVector() : data_( nullptr ), size_( 0 ), externally_allocated_( false ) {}

    ATLAS_HOST_DEVICE
    SVector( const T* data, const idx_t size ) : data_( data ), size_( size ), externally_allocated_( true ) {}

    ATLAS_HOST_DEVICE
    SVector( SVector const& other ) : data_( other.data_ ), size_( other.size_ ), externally_allocated_( true ) {}

    ATLAS_HOST_DEVICE
    SVector( SVector&& other ) :
        data_( other.data_ ),
        size_( other.size_ ),
        externally_allocated_( other.externally_allocated_ ) {}

    ATLAS_HOST_DEVICE
    SVector& operator=( SVector const& other ) {
        data_                 = other.data_;
        size_                 = other.size_;
        externally_allocated_ = true;
        return *this;
    }

    SVector& operator=( SVector&& other ) = default;

    ATLAS_HOST_DEVICE
    SVector( T* data, idx_t size ) : data_( data ), size_( size ), externally_allocated_( true ) {}

    SVector( idx_t N ) : data_( nullptr ), size_( N ), externally_allocated_( false ) {
        util::allocate_managedmem( data_, N );
    }
    ATLAS_HOST_DEVICE
    ~SVector() {
#ifndef __CUDA_ARCH__
        if ( !externally_allocated_ )
            util::delete_managedmem( data_ );
#endif
    }

    void insert( idx_t pos, idx_t dimsize ) {
        T* data;
        util::allocate_managedmem( data, size_ + dimsize );

        for ( idx_t c = 0; c < pos; ++c ) {
            data[c] = data_[c];
        }
        for ( idx_t c = pos; c < size_; ++c ) {
            data[c + dimsize] = data_[c];
        }

        T* oldptr = data_;
        data_     = data;
        util::delete_managedmem( oldptr );
        size_ += dimsize;
    }

    size_t footprint() const { return sizeof( T ) * size_; }

    ATLAS_HOST_DEVICE
    T* data() { return data_; }

    ATLAS_HOST_DEVICE
    T const* data() const { return data_; }

    ATLAS_HOST_DEVICE
    T& operator()( const idx_t idx ) {
        assert( data_ && idx < size_ );
        return data_[idx];
    }
    ATLAS_HOST_DEVICE
    T const& operator()( const idx_t idx ) const {
        assert( data_ && idx < size_ );
        return data_[idx];
    }

    ATLAS_HOST_DEVICE
    T& operator[]( const idx_t idx ) {
        assert( data_ && idx < size_ );
        return data_[idx];
    }
    ATLAS_HOST_DEVICE
    T const& operator[]( const idx_t idx ) const {
        assert( data_ && idx < size_ );
        return data_[idx];
    }

    ATLAS_HOST_DEVICE
    idx_t size() const { return size_; }

    void resize_impl( idx_t N ) {
        if ( N == size_ )
            return;

        T* d_ = nullptr;
        util::allocate_managedmem( d_, N );
        for ( idx_t c = 0; c < std::min( size_, N ); ++c ) {
            d_[c] = data_[c];
        }
        util::delete_managedmem( data_ );
        data_ = d_;
    }

    void resize( idx_t N ) {
        resize_impl( N );
        size_ = N;
    }


    void resize( idx_t N, T&& val ) {
        const int oldsize = size_;
        resize( N );
        for ( idx_t c = oldsize; c < size_; ++c ) {
            data_[c] = val;
        }
    }

private:
    T* data_;
    idx_t size_;
    bool externally_allocated_;
};

//------------------------------------------------------------------------------

}  // namespace array
}  // namespace atlas
