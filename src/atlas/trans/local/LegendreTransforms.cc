/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/trans/local/LegendreTransforms.h"
#include <cstddef>

namespace atlas {
namespace trans {

//-----------------------------------------------------------------------------

void invtrans_legendre( const size_t trc,    // truncation (in)
                        const size_t trcFT,  // truncation for Fourier transformation (in)
                        const size_t trcLP,  // truncation of Legendre polynomials data legpol. Needs to be >= trc (in)
                        const double legpol[],  // values of associated Legendre functions, size (trc+1)*trc/2 (in)
                        const int nb_fields,    // number of fields
                        const double spec[],    // spectral data, size (trc+1)*trc (in)
                        double leg_real[],      // values of associated Legendre functions, size (trc+1)*trc/2 (out)
                        double leg_imag[] )     // values of associated Legendre functions, size (trc+1)*trc/2 (out)
{
    // Legendre transformation:
    int k = 0, klp = 0;
    for ( int jm = 0; jm <= trcFT; ++jm ) {
        for ( int jfld = 0; jfld < nb_fields; ++jfld ) {
            leg_real[jm * nb_fields + jfld] = 0.;
            leg_imag[jm * nb_fields + jfld] = 0.;
        }
        for ( int jn = jm; jn <= trcLP; ++jn, ++klp ) {
            if ( jn <= trc ) {
                for ( int jfld = 0; jfld < nb_fields; ++jfld ) {
                    // not completely sure where this factor 2 comes from. One possible
                    // explanation:
                    // normalization of trigonometric functions in the spherical harmonics
                    // integral over square of trig function is 1 for m=0 and 0.5 (?) for
                    // m>0
                    leg_real[jm * nb_fields + jfld] += 2. * spec[( 2 * k ) * nb_fields + jfld] * legpol[klp];
                    leg_imag[jm * nb_fields + jfld] += 2. * spec[( 2 * k + 1 ) * nb_fields + jfld] * legpol[klp];
                }
                ++k;
            }
        }
    }
    // Undo factor 2 for (jm == 0)
    for ( int jfld = 0; jfld < nb_fields; ++jfld ) {
        leg_real[jfld] /= 2.;
        leg_imag[jfld] /= 2.;
    }
}

// --------------------------------------------------------------------------------------------------------------------

}  // namespace trans
}  // namespace atlas
