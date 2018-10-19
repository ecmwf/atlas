/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#pragma once

#include <vector>
#include "atlas/library/config.h"
#include "atlas/util/Config.h"

namespace atlas {

//---------------------------------------------------------------------------------------------------------------------

class Vertical {
public:
    template <typename vector_t>  // expect "vector_t::size()" and "vector_t::operator[]"
    Vertical( idx_t levels, const vector_t& z, const util::Config& config = util::NoConfig() );

    Vertical( const util::Config& config = util::NoConfig() );

public:
    idx_t k_begin() const { return k_begin_; }
    idx_t k_end() const { return k_end_; }
    bool boundaries() const { return boundaries_; }
    idx_t size() const { return size_; }

    template <typename Int>
    double operator()( const Int k ) const {
        return z_[k];
    }

    template <typename Int>
    double operator[]( const Int k ) const {
        return z_[k];
    }

private:
    bool boundaries_;
    idx_t k_begin_;
    idx_t k_end_;
    idx_t size_;
    std::vector<double> z_;
};

//---------------------------------------------------------------------------------------------------------------------

template <typename vector_t>
Vertical::Vertical( idx_t levels, const vector_t& z, const util::Config& config ) {
    size_       = levels;
    boundaries_ = config.getBool( "boundaries", false );
    k_begin_    = 0;
    k_end_      = size_;
    if ( boundaries_ ) {
        size_ += 2;
        ++k_begin_;
        k_end_ = size_ - 1;
    }
    ASSERT( size_ == static_cast<idx_t>( z.size() ) );
    z_.resize( size_ );
    for ( idx_t k = 0; k < size_; ++k ) {
        z_[k] = z[k];
    }
}

//---------------------------------------------------------------------------------------------------------------------

}  // namespace atlas
