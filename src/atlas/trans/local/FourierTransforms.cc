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

#include "atlas/trans/local/FourierTransforms.h"
#include <algorithm>
#include <cmath>
#include <iostream>

namespace atlas {
namespace trans {

//-----------------------------------------------------------------------------

void invtrans_fourier( const size_t trcFT,
                       const double lon,         // longitude in radians (in)
                       const int nb_fields,      // Number of fields
                       const double rlegReal[],  // values of associated Legendre
                                                 // functions, size (trc+1)*trc/2
                                                 // (in)
                       const double rlegImag[],  // values of associated Legendre
                                                 // functions, size (trc+1)*trc/2
                                                 // (in)
                       double rgp[] )            // gridpoint
{
    for ( int jfld = 0; jfld < nb_fields; ++jfld ) {
        rgp[jfld] = 0.;
    }
    // local Fourier transformation:
    // (FFT would be slower when computing the Fourier transformation for a single
    // point)
    for ( int jm = 0; jm <= trcFT; ++jm ) {
        const double cos = std::cos( jm * lon );
        const double sin = std::sin( jm * lon );
        for ( int jfld = 0; jfld < nb_fields; ++jfld ) {
            rgp[jfld] += cos * rlegReal[jm * nb_fields + jfld] - sin * rlegImag[jm * nb_fields + jfld];
        }
    }
}

int fourier_truncation( const int truncation, const int nx, const int nxmax, const int ndgl, const double lat,
                        const bool fullgrid ) {
    int linear_truncation = ndgl - 1;
    int trc               = truncation;
    if ( truncation >= linear_truncation || fullgrid ) {
        // linear
        trc = ( nx - 1 ) / 2;
    }
    else if ( truncation >= ndgl * 2 / 3 - 1 ) {
        // quadratic
        // trc = (nx-1)/(2+std::pow(std::cos(lat),2));
        trc = ( nx - 1 ) / ( 2 + 3 * ( linear_truncation - truncation ) / ndgl * std::pow( std::cos( lat ), 2 ) );
    }
    else {
        // cubic
        trc = ( nx - 1 ) / ( 2 + std::pow( std::cos( lat ), 2 ) ) - 1;
    }
    trc = std::min( truncation, trc );
    // std::cout << "truncation=" <<  truncation << " trc=" << trc << "
    // ndgl*2/3-1=" << ndgl*2/3-1 <<  " nx=" << nx << " nxmax=" << nxmax << "
    // latsin2=" << std::pow(std::cos(lat),2) << std::endl;
    return trc;
}

// --------------------------------------------------------------------------------------------------------------------

}  // namespace trans
}  // namespace atlas
