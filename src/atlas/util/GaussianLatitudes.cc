/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

/// @file   GaussianLatitudes.cc
/// @author Willem Deconinck
/// @date   Jan 2014

#include "atlas/util/GaussianLatitudes.h"
#include "atlas/grid/detail/spacing/gaussian/Latitudes.h"

namespace atlas {
namespace util {

void gaussian_latitudes_npole_equator( const size_t N, double latitudes[] ) {
    grid::spacing::gaussian::gaussian_latitudes_npole_equator( N, latitudes );
}

void gaussian_quadrature_npole_equator( const size_t N, double latitudes[], double weights[] ) {
    grid::spacing::gaussian::gaussian_quadrature_npole_equator( N, latitudes, weights );
}

void gaussian_latitudes_npole_spole( const size_t N, double latitudes[] ) {
    grid::spacing::gaussian::gaussian_latitudes_npole_spole( N, latitudes );
}

void gaussian_quadrature_npole_spole( const size_t N, double latitudes[], double weights[] ) {
    grid::spacing::gaussian::gaussian_quadrature_npole_spole( N, latitudes, weights );
}

}  // namespace util
}  // namespace atlas
