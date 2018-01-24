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

#include <cstddef>

namespace atlas {
namespace trans {

//-----------------------------------------------------------------------------
// Routine to compute the Legendre transformation
//
// Author:
// Andreas Mueller *ECMWF*
//
void invtrans_legendre(
        const size_t trc,       // truncation (in)
        const size_t trcFT,     // truncation for Fourier transformation (in)
        const double legpol[],  // values of associated Legendre functions, size (trc+1)*trc/2 (in)
        const int nb_fields,    // Number of fields
        const double spec[],    // spectral data, size (trc+1)*trc (in)
        double leg_real[],      // values of associated Legendre functions, size (trc+1)*trc/2 (out)
        double leg_imag[] );    // values of associated Legendre functions, size (trc+1)*trc/2 (out)

// --------------------------------------------------------------------------------------------------------------------

} // namespace trans
} // namespace atlas
