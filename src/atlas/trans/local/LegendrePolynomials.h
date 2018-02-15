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
// Routine to compute the Legendre polynomials in serial according to Belousov
// (using correction by Swarztrauber)
//
// Reference:
// S.L. Belousov, Tables of normalized associated Legendre Polynomials, Pergamon
// Press (1962)
// P.N. Swarztrauber, On computing the points and weights for Gauss-Legendre
// quadrature,
//      SIAM J. Sci. Comput. Vol. 24 (3) pp. 945-954 (2002)
//
// Author of Fortran version:
// Mats Hamrud, Philippe Courtier, Nils Wedi *ECMWF*
//
// Ported to C++ by:
// Andreas Mueller *ECMWF*
//
void compute_legendre_polynomials( const size_t trc,   // truncation (in)
                                   const double lat,   // latitude in radians (in)
                                   double legpol[] );  // values of associated
                                                       // Legendre functions, size
                                                       // (trc+1)*trc/2 (out)

// --------------------------------------------------------------------------------------------------------------------

}  // namespace trans
}  // namespace atlas
