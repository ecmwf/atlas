/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/grid/Vertical.h"

namespace atlas {

//---------------------------------------------------------------------------------------------------------------------

namespace {

std::vector<double> linspace( double start, double end, idx_t N, bool endpoint ) {
    std::vector<double> x_;
    x_.resize( N );

    double step;
    if ( endpoint && N > 1 )
        step = ( end - start ) / double( N - 1 );
    else
        step = ( end - start ) / double( N );

    for ( idx_t i = 0; i < N; ++i ) {
        x_[i] = start + i * step;
    }
    return x_;
}

bool get_boundaries( const util::Config& config ) {
    return config.getBool( "boundaries", false );
}

idx_t get_levels( const util::Config& config ) {
    return config.getInt( "levels", 0 );
}

idx_t get_size( const util::Config& config ) {
    idx_t levels = get_levels( config );
    idx_t size   = levels ? levels + 2 * idx_t{get_boundaries( config )} : 0;
    return size;
}

}  // namespace

//---------------------------------------------------------------------------------------------------------------------

Vertical::Vertical( const util::Config& config ) :
    Vertical( get_levels( config ), linspace( 0., 1., get_size( config ), true ), config ) {}

//---------------------------------------------------------------------------------------------------------------------

}  // namespace atlas
