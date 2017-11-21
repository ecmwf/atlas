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

double invtrans_fourier(
        const size_t trcFT,
        const double rlegReal[],  // values of associated Legendre functions, size (trc+1)*trc/2 (out)
        const double rlegImag[],  // values of associated Legendre functions, size (trc+1)*trc/2 (out)
        const double lon )        // longitude in radians (in)
{
    double result(0);
    // local Fourier transformation:
    // (FFT would be slower when computing the Fourier transformation for a single point)
    for( int jm=0; jm<=trcFT; ++jm ) {
        result += std::cos(jm*lon) * rlegReal[jm] - std::sin(jm*lon) * rlegImag[jm];
    }
    return result;
}

// --------------------------------------------------------------------------------------------------------------------

} // namespace trans
} // namespace atlas
