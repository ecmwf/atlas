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

#include <cmath>
#include <limits>

#include "atlas/array.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/trans/localopt/LegendrePolynomialsopt.h"

namespace atlas {
namespace trans {

//-----------------------------------------------------------------------------

void compute_legendre_polynomialsopt(
    const size_t trc,          // truncation (in)
    const int nlats,           // number of latitudes
    const double lats[],       // latitudes in radians (in)
    double leg_sym[],          // values of associated Legendre functions, symmetric part
    double leg_asym[],         // values of associated Legendre functions, asymmetric part
    size_t leg_start_sym[],    // start indices for different zonal wave numbers, symmetric part
    size_t leg_start_asym[] )  // start indices for different zonal wave numbers, asymmetric part
{
    auto legendre_size = [&]( int truncation ) { return ( truncation + 2 ) * ( truncation + 1 ) / 2; };
    array::ArrayT<double> zfn_( trc + 1, trc + 1 );
    array::ArrayView<double, 2> zfn = array::make_view<double, 2>( zfn_ );
    std::vector<double> legpol( legendre_size( trc ) );
    auto idxmn = [&]( int jm, int jn ) { return ( 2 * trc + 3 - jm ) * jm / 2 + jn - jm; };
    int iodd;

    // Compute coefficients for Taylor series in Belousov (19) and (21)
    // Belousov, Swarztrauber use zfn(0,0)=std::sqrt(2.)
    // IFS normalisation chosen to be 0.5*Integral(Pnm**2) = 1
    zfn( 0, 0 ) = 2.;
    for ( int jn = 1; jn <= trc; ++jn ) {
        double zfnn = zfn( 0, 0 );
        for ( int jgl = 1; jgl <= jn; ++jgl ) {
            zfnn *= std::sqrt( 1. - 0.25 / ( jgl * jgl ) );
        }
        iodd          = jn % 2;
        zfn( jn, jn ) = zfnn;
        for ( int jgl = 2; jgl <= jn - iodd; jgl += 2 ) {
            double zfjn = ( ( jgl - 1. ) * ( 2. * jn - jgl + 2. ) );  // new factor numerator
            double zfjd = ( jgl * ( 2. * jn - jgl + 1. ) );           // new factor denominator

            zfn( jn, jn - jgl ) = zfn( jn, jn - jgl + 2 ) * zfjn / zfjd;
        }
    }

    // Loop over latitudes:
    for ( int jlat = 0; jlat < nlats; ++jlat ) {
        {
            //ATLAS_TRACE( "compute Legendre polynomials" );
            // --------------------
            // 1. First two columns
            // --------------------
            double lat     = lats[jlat];
            double zdlx1   = ( M_PI_2 - lat );               // theta
            double zdlx    = std::cos( zdlx1 );              // cos(theta)
            double zdlsita = std::sqrt( 1. - zdlx * zdlx );  // sin(theta) (this is how trans library does it)

            legpol[idxmn( 0, 0 )] = 1.;

            double zdl1sita = 0.;
            // if we are less than 1 meter from the pole,
            if ( std::abs( zdlsita ) <= std::sqrt( std::numeric_limits<double>::epsilon() ) ) {
                zdlx    = 1.;
                zdlsita = 0.;
            }
            else {
                zdl1sita = 1. / zdlsita;
            }

            // ordinary Legendre polynomials from series expansion
            // ---------------------------------------------------

            // even N
            for ( int jn = 2; jn <= trc; jn += 2 ) {
                double zdlk   = 0.5 * zfn( jn, 0 );
                double zdlldn = 0.0;
                double zdsq   = 1. / std::sqrt( jn * ( jn + 1. ) );
                // represented by only even k
                for ( int jk = 2; jk <= jn; jk += 2 ) {
                    // normalised ordinary Legendre polynomial == \overbar{P_n}^0
                    zdlk = zdlk + zfn( jn, jk ) * std::cos( jk * zdlx1 );
                    // normalised associated Legendre polynomial == \overbar{P_n}^1
                    zdlldn = zdlldn + zdsq * zfn( jn, jk ) * jk * std::sin( jk * zdlx1 );
                }
                legpol[idxmn( 0, jn )] = zdlk;
                legpol[idxmn( 1, jn )] = zdlldn;
            }

            // odd N
            for ( int jn = 1; jn <= trc; jn += 2 ) {
                zfn( jn, 0 )  = 0.;
                double zdlk   = 0.;
                double zdlldn = 0.0;
                double zdsq   = 1. / std::sqrt( jn * ( jn + 1. ) );
                // represented by only even k
                for ( int jk = 1; jk <= jn; jk += 2 ) {
                    // normalised ordinary Legendre polynomial == \overbar{P_n}^0
                    zdlk = zdlk + zfn( jn, jk ) * std::cos( jk * zdlx1 );
                    // normalised associated Legendre polynomial == \overbar{P_n}^1
                    zdlldn = zdlldn + zdsq * zfn( jn, jk ) * jk * std::sin( jk * zdlx1 );
                }
                legpol[idxmn( 0, jn )] = zdlk;
                legpol[idxmn( 1, jn )] = zdlldn;
            }

            // --------------------------------------------------------------
            // 2. Diagonal (the terms 0,0 and 1,1 have already been computed)
            //    Belousov, equation (23)
            // --------------------------------------------------------------

            double zdls = zdl1sita * std::numeric_limits<double>::min();
            for ( int jn = 2; jn <= trc; ++jn ) {
                double sq = std::sqrt( ( 2. * jn + 1. ) / ( 2. * jn ) );

                legpol[idxmn( jn, jn )] = legpol[idxmn( jn - 1, jn - 1 )] * zdlsita * sq;
                if ( std::abs( legpol[idxmn( jn, jn )] ) < zdls ) legpol[idxmn( jn, jn )] = 0.0;
            }

            // ---------------------------------------------
            // 3. General recurrence (Belousov, equation 17)
            // ---------------------------------------------

            for ( int jn = 3; jn <= trc; ++jn ) {
                for ( int jm = 2; jm < jn; ++jm ) {
                    double cn =
                        ( ( 2. * jn + 1. ) * ( jn + jm - 3. ) * ( jn + jm - 1. ) );     // numerator of c in Belousov
                    double cd = ( ( 2. * jn - 3. ) * ( jn + jm - 2. ) * ( jn + jm ) );  // denominator of c in Belousov
                    double dn =
                        ( ( 2. * jn + 1. ) * ( jn - jm + 1. ) * ( jn + jm - 1. ) );     // numerator of d in Belousov
                    double dd = ( ( 2. * jn - 1. ) * ( jn + jm - 2. ) * ( jn + jm ) );  // denominator of d in Belousov
                    double en = ( ( 2. * jn + 1. ) * ( jn - jm ) );                     // numerator of e in Belousov
                    double ed = ( ( 2. * jn - 1. ) * ( jn + jm ) );                     // denominator of e in Belousov

                    legpol[idxmn( jm, jn )] = std::sqrt( cn / cd ) * legpol[idxmn( jm - 2, jn - 2 )] -
                                              std::sqrt( dn / dd ) * legpol[idxmn( jm - 2, jn - 1 )] * zdlx +
                                              std::sqrt( en / ed ) * legpol[idxmn( jm, jn - 1 )] * zdlx;
                }
            }
        }

        {
            //ATLAS_TRACE( "add to global arrays" );

            // take factor 2 for m > 0 into account:
            for ( int jm = 1; jm <= trc; ++jm ) {
                for ( int jn = jm; jn <= trc; ++jn ) {
                    //legpol[idxmn( jm, jn )] *= 2.;
                }
            }
            for ( int jm = 0; jm <= trc; jm++ ) {
                int is1 = 0, ia1 = 0;
                for ( int jn = jm; jn <= trc; jn++ ) {
                    if ( ( jn - jm ) % 2 == 0 ) { is1++; }
                    else {
                        ia1++;
                    }
                }
                int is2 = 0, ia2 = 0;
                // the choice between the following two code lines determines whether
                // total wavenumbers are summed in an ascending or descending order.
                // The trans library in IFS uses descending order because it should
                // be more accurate (higher wavenumbers have smaller contributions).
                // This also needs to be changed when splitting the spectral data in
                // TransLocalopt::invtrans_uv!
                //for ( int jn = jm; jn <= trc; jn++ ) {
                for ( int jn = trc; jn >= jm; jn-- ) {
                    if ( ( jn - jm ) % 2 == 0 ) {
                        int is      = leg_start_sym[jm] + jlat + nlats * is2++;
                        leg_sym[is] = legpol[idxmn( jm, jn )];
                    }
                    else {
                        int ia       = leg_start_asym[jm] + jlat + nlats * ia2++;
                        leg_asym[ia] = legpol[idxmn( jm, jn )];
                    }
                }
            }
        }
    }
}

// --------------------------------------------------------------------------------------------------------------------

}  // namespace trans
}  // namespace atlas
