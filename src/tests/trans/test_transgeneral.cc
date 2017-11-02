/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */


#include <algorithm>

#include "atlas/library/Library.h"
#include "atlas/field/FieldSet.h"
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/functionspace/Spectral.h"
#include "atlas/functionspace/StructuredColumns.h"
#include "atlas/grid/Distribution.h"
#include "atlas/grid/Partitioner.h"
#include "atlas/grid.h"
#include "atlas/grid/detail/partitioner/EqualRegionsPartitioner.h"
#include "atlas/grid/detail/partitioner/TransPartitioner.h"
#include "atlas/meshgenerator/StructuredMeshGenerator.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/output/Gmsh.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/trans/Trans.h"
#include "atlas/array/MakeView.h"
#include "transi/trans.h"

#include "tests/AtlasTestEnvironment.h"
#include "eckit/testing/Test.h"

#include <iomanip>

using namespace eckit;
using namespace eckit::testing;

using atlas::array::Array;
using atlas::array::ArrayView;
using atlas::array::make_view;

namespace atlas {
namespace test {

//-----------------------------------------------------------------------------

struct AtlasTransEnvironment : public AtlasTestEnvironment {
       AtlasTransEnvironment(int argc, char * argv[]) : AtlasTestEnvironment(argc, argv) {
         if( parallel::mpi::comm().size() == 1 )
           trans_use_mpi(false);
         trans_init();
       }

      ~AtlasTransEnvironment() {
         trans_finalize();
       }
};

//-----------------------------------------------------------------------------
// Routine to compute the Legendre polynomials in serial according to Belousov
// (using correction by Swarztrauber)
//
// Reference:
// S.L. Belousov, Tables of normalized associated Legendre Polynomials, Pergamon Press (1962)
// P.N. Swarztrauber, On computing the points and weights for Gauss-Legendre quadrature,
//      SIAM J. Sci. Comput. Vol. 24 (3) pp. 945-954 (2002)
//
// Author of Fortran version:
// Mats Hamrud, Philippe Courtier, Nils Wedi *ECMWF*
//
// Ported to C++ by:
// Andreas Mueller *ECMWF*
//
void compute_legendre(
        const size_t knsmax,               // truncation + 1 (in)
        double& lat,                       // latitude in degree (in)
        array::ArrayView<double,1>& zlfpol)// values of associated Legendre functions, size (knsmax+1)*knsmax/2 (out)
{
    std::ostream& out = Log::info(); // just for debugging
    array::ArrayT<int> idxmn_(knsmax+1,knsmax+1);
    array::ArrayView<int,2> idxmn = make_view<int,2>(idxmn_);
    int j = 0;
    for( int jm=0; jm<=knsmax; ++jm ) {
        for( int jn=jm; jn<=knsmax; ++jn ) {
            idxmn(jm,jn) = j++;
        }
    }

    array::ArrayT<double> zfn_(knsmax+1,knsmax+1);
    array::ArrayView<double,2> zfn = make_view<double,2>(zfn_);

    int iodd;

    // Belousov, Swarztrauber use zfn(0,0)=std::sqrt(2.)
    // IFS normalisation chosen to be 0.5*Integral(Pnm**2) = 1
    zfn(0,0) = 2.;
    for( int jn=1; jn<=knsmax; ++jn ) {
        double zfnn = zfn(0,0);
        for( int jgl=1; jgl<=jn; ++jgl) {
            zfnn *= std::sqrt(1.-0.25/(jgl*jgl));
        }
        iodd = jn % 2;
        zfn(jn,jn)=zfnn;
        for( int jgl=2; jgl<=jn-iodd; jgl+=2 ) {
            zfn(jn,jn-jgl) = zfn(jn,jn-jgl+2) * ((jgl-1.)*(2.*jn-jgl+2.)) / (jgl*(2.*jn-jgl+1.));
        }
    }

    zlfpol.assign(0.); // just for debugging (shouldn't be necessary because not in trans library)

    // --------------------
    // 1. First two columns
    // --------------------
    double zdlx1 = (90.-lat)*util::Constants::degreesToRadians(); // theta
    double zdlx = std::cos(zdlx1); // cos(theta)
    double zdlsita = std::sqrt(1.-zdlx*zdlx); // sin(theta) (this is how trans library does it)

    zlfpol(0,0) = 1.;
    double zdl1sita = 0.;

    // if we are less than 1 meter from the pole,
    if( std::abs(zdlsita) <= std::sqrt(std::numeric_limits<double>::epsilon()) )
    {
        zdlx = 1.;
        zdlsita = 0.;
    } else {
        zdl1sita = 1./zdlsita;
    }

    // ordinary Legendre polynomials from series expansion
    // ---------------------------------------------------

    // even N
    for( int jn=2; jn<=knsmax; jn+=2 ) {
        double zdlk = 0.5*zfn(jn,0);
        double zdlldn = 0.0;
        // represented by only even k
        for (int jk=2; jk<=jn; jk+=2 ) {
            // normalised ordinary Legendre polynomial == \overbar{P_n}^0
            zdlk = zdlk + zfn(jn,jk)*std::cos(jk*zdlx1);
            // normalised associated Legendre polynomial == \overbar{P_n}^1
            zdlldn = zdlldn + 1./std::sqrt(jn*(jn+1.))*zfn(jn,jk)*jk*std::sin(jk*zdlx1);
        }
        zlfpol(idxmn(0,jn)) = zdlk;
        zlfpol(idxmn(1,jn)) = zdlldn;
    }

    // odd N
    for( int jn=1; jn<=knsmax; jn+=2 ) {
        double zdlk = 0.5*zfn(jn,0);
        double zdlldn = 0.0;
        // represented by only even k
        for (int jk=1; jk<=jn; jk+=2 ) {
            // normalised ordinary Legendre polynomial == \overbar{P_n}^0
            zdlk = zdlk + zfn(jn,jk)*std::cos(jk*zdlx1);
            // normalised associated Legendre polynomial == \overbar{P_n}^1
            zdlldn = zdlldn + 1./std::sqrt(jn*(jn+1.))*zfn(jn,jk)*jk*std::sin(jk*zdlx1);
        }
        zlfpol(idxmn(0,jn)) = zdlk;
        zlfpol(idxmn(1,jn)) = zdlldn;
    }

    // --------------------------------------------------------------
    // 2. Diagonal (the terms 0,0 and 1,1 have already been computed)
    //    Belousov, equation (23)
    // --------------------------------------------------------------

    double zdls = zdl1sita*std::numeric_limits<double>::min();
    for( int jn=2; jn<=knsmax; ++jn ) {
        zlfpol(idxmn(jn,jn)) = zlfpol(idxmn(jn-1,jn-1))*zdlsita*std::sqrt((2.*jn+1.)/(2.*jn));
        if( std::abs(zlfpol(idxmn(jn,jn))) < zdls ) zlfpol(idxmn(jn,jn)) = 0.0;
    }

    // ---------------------------------------------
    // 3. General recurrence (Belousov, equation 17)
    // ---------------------------------------------

    for( int jn=3; jn<=knsmax; ++jn ) {
        for( int jm=2; jm<jn; ++jm ) {
            zlfpol(idxmn(jm,jn)) =
                    std::sqrt(((2.*jn+1.)*(jn+jm-1.)*(jn+jm-3.))/
                              ((2.*jn-3.)*(jn+jm   )*(jn+jm-2.))) * zlfpol(idxmn(jm-2,jn-2))
                  - std::sqrt(((2.*jn+1.)*(jn+jm-1.)*(jn-jm+1.))/
                              ((2.*jn-1.)*(jn+jm   )*(jn+jm-2.))) * zlfpol(idxmn(jm-2,jn-1)) * zdlx
                  + std::sqrt(((2.*jn+1.)*(jn-jm   )           )/
                              ((2.*jn-1.)*(jn+jm   )           )) * zlfpol(idxmn(jm,jn-1  )) * zdlx;
        }
    }

}

//-----------------------------------------------------------------------------
// Routine to compute the spectral transform by using a local Fourier transformation
// for a single point
//
// Author:
// Andreas Mueller *ECMWF*
//
double spectral_transform_point(
        const size_t knsmax,               // truncation + 1 (in)
        const size_t knsmaxFT,             // truncation + 1 for Fourier transformation (in)
        double& lon,                       // longitude in degree (in)
        array::ArrayView<double,1>& zlfpol,// values of associated Legendre functions at desired latitude, size (knsmax+1)*knsmax/2 (in)
        array::ArrayView<double,1>& rspecg)// spectral data, size (knsmax+1)*knsmax (in)
{
    double result = 0., zfac = 1.;
    int k = 0;
    double lonRad = lon * util::Constants::degreesToRadians();
    for( int jm=0; jm<=knsmax; ++jm ) {
        // Legendre transformation:
        double rlegReal = 0., rlegImag = 0.;
        for( int jn=jm; jn<=knsmax; ++jn ) {
            if( jm>0 ) zfac = 2.;
            // not completely sure where this zfac comes from. One possible explanation:
            // normalization of trigonometric functions in the spherical harmonics
            // integral over square of trig function is 1 for m=0 and 0.5 (?) for m>0
            rlegReal += rspecg(2*k)   * zfac * zlfpol(k);
            rlegImag += rspecg(2*k+1) * zfac * zlfpol(k);
            //Log::info() << zfac * zlfpol(k) << std::endl; // just for debugging
            k++;
        }
        // local Fourier transformation:
        // (FFT would be slower when computing the Fourier transformation for a single point)
        if( jm<=knsmaxFT ) result += std::cos(jm*lonRad) * rlegReal - std::sin(jm*lonRad) * rlegImag;
    }
    return result;
}


//-----------------------------------------------------------------------------

CASE( "test_transgeneral_legendrepolynomials" )
{
  std::ostream& out = Log::info(); // just for debugging
  out << "test_transgeneral_legendrepolynomials" << std::endl;
/*
  Grid g( "O10" );
  trans::Trans trans(g,1279);
*/
  int trp1 = 1280; // truncation + 1
  int N = (trp1+2)*(trp1+1)/2;
  atlas::array::ArrayT<double> zlfpol_(N);
  atlas::array::ArrayView<double,1> zlfpol = make_view<double,1>(zlfpol_);

  double lat = 50.;
  lat = std::acos(0.99312859918509488)*util::Constants::radiansToDegrees();
  compute_legendre(trp1, lat, zlfpol);
/*
  out << "zlfpol after compute legendre:" << std::endl;
  int idx, col = 8;
  for( int jn=0; jn<std::ceil(1.*N/col); ++jn ) {
      for( int jm=0; jm<col; ++jm ) {
          idx = col*jn+jm;
          if( idx < N ) out << std::setw(15) << std::left << zlfpol(idx);
      }
      out << std::endl;
  }
*/

  //EXPECT(eckit::testing::make_view(arrv.data(),arrv.data()+2*f.N) == eckit::testing::make_view(arr_c,arr_c+2*f.N)); break; }

}

//-----------------------------------------------------------------------------

CASE( "test_transgeneral_pointtrans" )
{
  std::ostream& out = Log::info();
  Log::info() << "test_transgeneral_pointtrans" << std::endl;

  int trp1 = 10; // truncation + 1
  int N = (trp1+2)*(trp1+1)/2;

  atlas::array::ArrayT<double> rspecg_(2*N);
  atlas::array::ArrayView<double,1> rspecg = make_view<double,1>(rspecg_);
  rspecg.assign( { // copy and paste from Python output from print repr(data) for geopotential of T1279 truncated to T10
                           2.27058862e+03,   0.00000000e+00,   2.64344482e+02,
                           0.00000000e+00,   1.04363721e+03,   0.00000000e+00,
                          -1.38711157e+03,   0.00000000e+00,   6.30294678e+02,
                           0.00000000e+00,  -1.47959766e+03,   0.00000000e+00,
                           1.25787305e+03,   0.00000000e+00,  -3.47238281e+02,
                           0.00000000e+00,   6.09284912e+02,   0.00000000e+00,
                          -5.79417480e+02,   0.00000000e+00,  -7.42720642e+01,
                           0.00000000e+00,   4.70171387e+02,  -4.99296387e+02,
                          -4.64239407e+00,  -4.20254883e+02,  -2.01318069e+02,
                          -4.17947510e+02,  -8.64668579e+01,   6.79094482e+02,
                           8.39252777e+01,   5.26367493e+01,  -7.47528839e+01,
                           7.56367920e+02,   2.95226318e+02,  -4.45547119e+02,
                           6.53360596e+01,   3.04475098e+02,   1.98545914e+02,
                          -6.05724854e+02,  -4.66925335e+00,   4.36788086e+02,
                          -4.38317627e+02,  -2.01735626e+02,  -6.73341553e+02,
                          -3.45433105e+02,  -7.00174316e+02,  -1.59601822e+01,
                          -1.14086212e+02,   1.66471054e+02,   2.38090469e+02,
                           1.47945251e+02,   5.53364014e+02,   1.67163818e+02,
                           8.92426300e+01,   1.93021362e+02,   3.87085419e+01,
                           7.25012970e+01,  -3.77425781e+02,   1.46001043e+01,
                           2.06437378e+01,  -2.54263626e+02,   2.88258545e+02,
                           4.34750977e+02,   2.13519592e+02,   3.96897217e+02,
                           8.74137115e+01,   7.21976471e+01,   1.45806274e+02,
                          -1.06001190e+02,   4.55372467e+01,  -1.79682510e+02,
                           4.84295959e+01,  -1.41918839e+02,  -1.50270279e+02,
                           2.25189957e+02,   1.10319427e+02,  -4.35088135e+02,
                           5.34815430e+02,   3.42563721e+02,   4.42061523e+02,
                           1.75658325e+02,   1.22033813e+02,  -2.49562073e+01,
                          -1.15247650e+02,   3.08883514e+01,  -3.12923828e+02,
                          -1.02068848e+02,  -3.29612549e+02,  -1.96804504e+02,
                          -1.12869827e+02,  -3.42539062e+02,  -2.32903732e+02,
                          -9.58003235e+01,  -7.35217896e+01,  -3.16965576e+02,
                          -1.24462708e+02,  -3.18577637e+02,  -1.14058228e+02,
                          -2.69070282e+01,  -3.63590851e+01,   6.86552277e+01,
                          -1.93415085e+02,  -3.96717262e+00,  -1.63823044e+02,
                           6.96005821e-01,  -6.39315796e+01,  -9.11142426e+01,
                          -1.09771667e+02,  -1.34256149e+02,  -1.35531940e+01,
                           1.38606615e+01,  -1.35011963e+02,   2.22399918e+02,
                           3.54877930e+02,   1.22672028e+02,   1.83927261e+02,
                           2.95129639e+02,   8.63545532e+01,   2.30139908e+02,
                          -1.14560532e+02,   6.74462585e+01,   3.10108154e+02,
                           9.13583469e+00,   1.77570038e+01,   1.12481117e+01,
                          -2.94228516e+01,  -2.62760925e+01,   7.95001831e+01,
                          -8.78986206e+01,   1.31246429e+02,   6.75210419e+01
               }
               );

  double lat = 27.9878, lon = 86.9250; // latitude and longitude in degree

  atlas::array::ArrayT<double> zlfpol_(N);
  atlas::array::ArrayView<double,1> zlfpol = make_view<double,1>(zlfpol_);

  // compute associated Legendre functions:
  compute_legendre(trp1, lat, zlfpol);

  // perform spectral transform:
  double result = spectral_transform_point(trp1, trp1, lon, zlfpol, rspecg);

  // output result:
  out << "result: " << result << std::endl;

}

//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas


int main(int argc, char **argv) {
    atlas::test::AtlasTransEnvironment env( argc, argv );
    return run_tests ( argc, argv, false );
}
