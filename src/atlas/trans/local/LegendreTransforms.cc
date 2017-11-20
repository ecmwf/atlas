/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <cstddef>
#include "atlas/trans/local/LegendreTransforms.h"

namespace atlas {
namespace trans {

//-----------------------------------------------------------------------------

void invtrans_legendre(
        const size_t trc,       // truncation (in)
        const size_t trcFT,     // truncation for Fourier transformation (in)
        const double legpol[],  // values of associated Legendre functions, size (trc+1)*trc/2 (in)
        const double spec[],    // spectral data, size (trc+1)*trc (in)
        double leg_real[],      // values of associated Legendre functions, size (trc+1)*trc/2 (out)
        double leg_imag[] )     // values of associated Legendre functions, size (trc+1)*trc/2 (out)
{
    // Legendre transformation:
    int k = 0;
    for( int jm=0; jm<=trcFT; ++jm ) {
        leg_real[jm] = 0.;
        leg_imag[jm] = 0.;
        for( int jn=jm; jn<=trc; ++jn, ++k ) {
            // not completely sure where this factor 2 comes from. One possible explanation:
            // normalization of trigonometric functions in the spherical harmonics
            // integral over square of trig function is 1 for m=0 and 0.5 (?) for m>0
            leg_real[jm] += 2. * spec[2*k]   * legpol[k];
            leg_imag[jm] += 2. * spec[2*k+1] * legpol[k];
        }
    }
    // Undo factor 2 for (jm == 0)
    leg_real[0] /= 2.;
    leg_imag[0] /= 2.;

}

// --------------------------------------------------------------------------------------------------------------------

} // namespace trans
} // namespace atlas
