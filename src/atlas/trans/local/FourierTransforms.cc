/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <cmath>
#include "atlas/trans/local/FourierTransforms.h"

namespace atlas {
namespace trans {

//-----------------------------------------------------------------------------

void invtrans_fourier(
        const size_t trcFT,
        const double lon,         // longitude in radians (in)
        const int nb_fields,      // Number of fields
        const double rlegReal[],  // values of associated Legendre functions, size (trc+1)*trc/2 (in)
        const double rlegImag[],  // values of associated Legendre functions, size (trc+1)*trc/2 (in)
        double rgp[] )            // gridpoint
{
    for( int jfld=0; jfld<nb_fields; ++jfld ) {
      rgp[jfld] = 0.;
    }
    // local Fourier transformation:
    // (FFT would be slower when computing the Fourier transformation for a single point)
    for( int jm=0; jm<=trcFT; ++jm ) {
        const double cos = std::cos(jm*lon);
        const double sin = std::sin(jm*lon);
        for( int jfld=0; jfld<nb_fields; ++jfld ) {
            rgp[jfld] += cos * rlegReal[jm*nb_fields+jfld] - sin * rlegImag[jm*nb_fields+jfld];
        }
    }
}

// --------------------------------------------------------------------------------------------------------------------

} // namespace trans
} // namespace atlas
