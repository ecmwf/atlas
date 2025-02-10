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
#include <vector>

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
void compute_zfn(const int trc, double zfn[]);


// Workspace to avoid repeated allocations
struct LegendrePolynomialsWorkspace {
    LegendrePolynomialsWorkspace(int trc) {
        vsin.reserve(trc+1);
        vcos.reserve(trc+1);
    }
    std::vector<double> vsin;
    std::vector<double> vcos;
};

void compute_legendre_polynomials_lat(const int trc,     // truncation (in)
                                      const double lat,  // latitude in radians (in)
                                      double legpol[],   // legendre polynomials
                                      double zfn[],
                                      LegendrePolynomialsWorkspace& w); // workspace to avoid allocations

void compute_legendre_polynomials(
    const int trc,             // truncation (in)
    const int nlats,           // number of latitudes
    const double lats[],       // latitudes in radians (in)
    double legendre_sym[],     // values of associated Legendre functions, symmetric part
    double legendre_asym[],    // values of associated Legendre functions, asymmetric part
    size_t leg_start_sym[],    // start indices for different zonal wave numbers, symmetric part
    size_t leg_start_asym[]);  // start indices for different zonal wave numbers, asymmetric part

void compute_legendre_polynomials_all(const int trc,        // truncation (in)
                                      const int nlats,      // number of latitudes
                                      const double lats[],  // latitudes in radians (in)
                                      double legendre[]);   // legendre polynomials for all latitudes

// --------------------------------------------------------------------------------------------------------------------

}  // namespace trans
}  // namespace atlas
