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

#include <cstddef>

namespace atlas {
namespace trans {

//-----------------------------------------------------------------------------
// Routine to compute the local Fourier transformation
//
// Author:
// Andreas Mueller *ECMWF*
//

void invtrans_fourier( const size_t trcFT,
                       const double lon,         // longitude in radians (in)
                       const int nb_fields,      // Number of fields
                       const double rlegReal[],  // values of associated Legendre
                                                 // functions, size (trc+1)*trc/2
                                                 // (in)
                       const double rlegImag[],  // values of associated Legendre
                                                 // functions, size (trc+1)*trc/2
                                                 // (in)
                       double rgp[] );           // gridpoint

int fourier_truncation( const int truncation, const int nx, const int nxmax, const int ndgl, const double lat,
                        const bool fullgrid );

// --------------------------------------------------------------------------------------------------------------------

}  // namespace trans
}  // namespace atlas
