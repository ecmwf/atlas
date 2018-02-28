/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor
 * does it submit to any jurisdiction.
 */

#include <algorithm>
#include <cmath>
#include <iostream>

#include "atlas/trans/localopt/FourierTransformsopt.h"

namespace atlas {
namespace trans {

//-----------------------------------------------------------------------------

void invtrans_fourieropt( const size_t trcFT,
                          const double lon,         // longitude in radians (in)
                          const int nb_fields,      // Number of fields
                          const double rlegReal[],  // associated Legendre functions, size (trc+1)*trc/2 (in)
                          const double rlegImag[],  // associated Legendre functions, size (trc+1)*trc/2 (in)
                          double rgp[] )            // gridpoint
{
    for ( int jfld = 0; jfld < nb_fields; ++jfld ) {
        rgp[jfld] = 0.;
    }
    // local Fourier transformation:
    for ( int jm = 0; jm <= trcFT; ++jm ) {
        const double cos = std::cos( jm * lon );
        const double sin = std::sin( jm * lon );
        for ( int jfld = 0; jfld < nb_fields; ++jfld ) {
            double real = cos * rlegReal[jm * nb_fields + jfld];
            double imag = sin * rlegImag[jm * nb_fields + jfld];
            rgp[jfld] += real - imag;
        }
    }
}

int fourier_truncationopt( const int truncation,    // truncation
                           const int nx,            // number of longitudes
                           const int nxmax,         // maximum nx
                           const int ndgl,          // number of latitudes
                           const double lat,        // latitude in radian
                           const bool fullgrid ) {  // regular grid
    int trc     = truncation;
    int trclin  = ndgl - 1;
    int trcquad = ndgl * 2 / 3 - 1;
    if ( truncation >= trclin || fullgrid ) {
        // linear
        trc = ( nx - 1 ) / 2;
    }
    else if ( truncation >= trcquad ) {
        // quadratic
        double weight = 3 * ( trclin - truncation ) / ndgl;
        double sqcos  = std::pow( std::cos( lat ), 2 );

        trc = ( nx - 1 ) / ( 2 + weight * sqcos );
    }
    else {
        // cubic
        double sqcos = std::pow( std::cos( lat ), 2 );

        trc = ( nx - 1 ) / ( 2 + sqcos ) - 1;
    }
    trc = std::min( truncation, trc );
    return trc;
}

// --------------------------------------------------------------------------------------------------------------------

}  // namespace trans
}  // namespace atlas
