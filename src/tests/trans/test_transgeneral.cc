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
#include "atlas/runtime/Trace.h"
#include "transi/trans.h"

#include "tests/AtlasTestEnvironment.h"
#include "eckit/testing/Test.h"

#include <iomanip>
#include <chrono>

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
// Routine to compute the Legendre transformation
//
// Author:
// Andreas Mueller *ECMWF*
//
void legendre_transform(
        const size_t knsmax,                 // truncation + 1 (in)
        const size_t knsmaxFT,               // truncation + 1 for Fourier transformation (in)
        array::ArrayView<double,1>& rlegReal,// values of associated Legendre functions, size (knsmax+1)*knsmax/2 (out)
        array::ArrayView<double,1>& rlegImag,// values of associated Legendre functions, size (knsmax+1)*knsmax/2 (out)
        array::ArrayView<double,1>& zlfpol,  // values of associated Legendre functions, size (knsmax+1)*knsmax/2 (out)
        double rspecg[])                     // spectral data, size (knsmax+1)*knsmax (in)
{
    // Legendre transformation:
    rlegReal.assign(0.); rlegImag.assign(0.);
    double zfac = 1.;
    int k = 0;
    for( int jm=0; jm<=knsmaxFT; ++jm ) {
        for( int jn=jm; jn<=knsmax; ++jn ) {
            if( jm>0 ) zfac = 2.;

            // not completely sure where this zfac comes from. One possible explanation:
            // normalization of trigonometric functions in the spherical harmonics
            // integral over square of trig function is 1 for m=0 and 0.5 (?) for m>0
            rlegReal[jm] += rspecg[2*k]   * zfac * zlfpol(k);
            rlegImag[jm] += rspecg[2*k+1] * zfac * zlfpol(k);
            //Log::info() << zfac * zlfpol(k) << std::endl; // just for debugging
            k++;
        }
    }
}

//-----------------------------------------------------------------------------
// Routine to compute the local Fourier transformation
//
// Author:
// Andreas Mueller *ECMWF*
//
double fourier_transform(
        const size_t knsmaxFT,
        array::ArrayView<double,1>& rlegReal,// values of associated Legendre functions, size (knsmax+1)*knsmax/2 (out)
        array::ArrayView<double,1>& rlegImag,// values of associated Legendre functions, size (knsmax+1)*knsmax/2 (out)
        double lon)
{
    // local Fourier transformation:
    // (FFT would be slower when computing the Fourier transformation for a single point)
    double lonRad = lon * util::Constants::degreesToRadians(), result = 0.;
    for( int jm=0; jm<=knsmaxFT; ++jm ) {
        result += std::cos(jm*lonRad) * rlegReal(jm) - std::sin(jm*lonRad) * rlegImag(jm);
    }
    return result;
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
        double lon,                        // longitude in degree (in)
        double lat,                        // latitude in degree (in)
        double rspecg[])                   // spectral data, size (knsmax+1)*knsmax (in)
{
    std::ostream& out = Log::info(); // just for debugging
    int N = (knsmax+2)*(knsmax+1)/2;
    ATLAS_TRACE();
    atlas::array::ArrayT<double> zlfpol_(N);
    atlas::array::ArrayView<double,1> zlfpol = make_view<double,1>(zlfpol_);

    atlas::array::ArrayT<double> rlegReal_(knsmaxFT+1);
    atlas::array::ArrayView<double,1> rlegReal = make_view<double,1>(rlegReal_);

    atlas::array::ArrayT<double> rlegImag_(knsmaxFT+1);
    atlas::array::ArrayView<double,1> rlegImag = make_view<double,1>(rlegImag_);

    // Legendre transform:
    compute_legendre(knsmax, lat, zlfpol);
    legendre_transform(knsmax, knsmaxFT, rlegReal, rlegImag, zlfpol, rspecg);

    // Fourier transform:
    return fourier_transform(knsmaxFT, rlegReal, rlegImag, lon);
}


//-----------------------------------------------------------------------------
// Routine to compute the spectral transform by using a local Fourier transformation
// for a grid (same latitude for all longitudes, allows to compute Legendre functions
// once for all longitudes)
//
// Author:
// Andreas Mueller *ECMWF*
//
void spectral_transform_grid(
        const size_t knsmax,  // truncation + 1 (in)
        const size_t knsmaxFT,// truncation + 1 for Fourier transformation (in)
        Grid grid,            // call with something like Grid("O32")
        double rspecg[],      // spectral data, size (knsmax+1)*knsmax (in)
        double rgp[],         // resulting grid point data (out)
        bool pointwise)       // use point function for unstructured mesh for testing purposes
{
    std::ostream& out = Log::info(); // just for debugging
    int N = (knsmax+2)*(knsmax+1)/2;
    ATLAS_TRACE();
    atlas::array::ArrayT<double> zlfpol_(N);
    atlas::array::ArrayView<double,1> zlfpol = make_view<double,1>(zlfpol_);

    atlas::array::ArrayT<double> rlegReal_(knsmaxFT+1);
    atlas::array::ArrayView<double,1> rlegReal = make_view<double,1>(rlegReal_);

    atlas::array::ArrayT<double> rlegImag_(knsmaxFT+1);
    atlas::array::ArrayView<double,1> rlegImag = make_view<double,1>(rlegImag_);

    for( int jm=0; jm<grid.size(); jm++) rgp[jm] = 0.;

    if( grid::StructuredGrid(grid) ) {
        grid::StructuredGrid g(grid);
        int idx = 0;
        for( size_t j=0; j<g.ny(); ++j ) {
            double lat = g.y(j);

            // Legendre transform:
            compute_legendre(knsmax, lat, zlfpol);
            legendre_transform(knsmax, knsmaxFT, rlegReal, rlegImag, zlfpol, rspecg);

            for( size_t i=0; i<g.nx(j); ++i ) {
                double lon = g.x(i,j);
                // Fourier transform:
                rgp[idx++] = fourier_transform(knsmaxFT, rlegReal, rlegImag, lon);
            }
        }
    } else {
        int idx = 0;
        for( PointXY p: grid.xy()) {
            double lon = p.x(), lat = p.y();
            if( pointwise ) {
                // alternative for testing: use spectral_transform_point function:
                rgp[idx++] = spectral_transform_point(knsmax, knsmaxFT, lon, lat, rspecg);
            } else {
                // Legendre transform:
                compute_legendre(knsmax, lat, zlfpol);
                legendre_transform(knsmax, knsmaxFT, rlegReal, rlegImag, zlfpol, rspecg);

                // Fourier transform:
                rgp[idx++] = fourier_transform(knsmaxFT, rlegReal, rlegImag, lon);
            }
        }
    }
}

//-----------------------------------------------------------------------------
// Routine to compute the spectral transform by computing spherical harmonics analytically
// (so far only wave number 3)
//
// Author:
// Andreas Mueller *ECMWF*
//
double sphericalharmonics_analytic(
        double n,             // total wave number (implemented so far for n<4
        double m,             // zonal wave number (implemented so far for m<4, m<n
        bool real,            // true: test real part, false: test imaginary part
        double lon,           // longitude in degree
        double lat)           // latitude in degree
{
    double lonRad = lon * util::Constants::degreesToRadians();
    double latRad = lat * util::Constants::degreesToRadians();
    double latsin = std::sin(latRad), latcos = std::cos(latRad);
    // Fourier part of the spherical harmonics:
    double rft = 1.; // not sure why I need a minus here
    if( m>0 ) rft *= 2.; // the famous factor 2 that noone really understands
    if( real ) {
        rft *= std::cos(m*lonRad);
    } else {
        rft *= -std::sin(m*lonRad);
    }
    // Legendre part of the spherical harmonics (following http://mathworld.wolfram.com/SphericalHarmonic.html
    // multiplied with -2*sqrt(pi) due to different normalization and different coordinates):
    // (can also be computed on http://www.wolframalpha.com with:
    // LegendreP[n, m, x]/Sqrt[1/2*Integrate[LegendreP[n, m, y]^2, {y, -1, 1}]])
    if ( m==0 && n==0 )
        return rft;
    if ( m==0 && n==1 )
        return std::sqrt(3.)*latsin*rft;
    if ( m==0 && n==2 )
        return std::sqrt(5.)/2.*(3.*latsin*latsin-1.)*rft; // shouldn't this be minus?
    if ( m==0 && n==3 )
        return std::sqrt(7.)/2.*(5.*latsin*latsin-3.)*latsin*rft; // shouldn't this be minus?
    if ( m==1 && n==1 )
        return std::sqrt(3./2.)*latcos*rft; // shouldn't this be minus?
    if ( m==1 && n==2 )
        return std::sqrt(15./2.)*latsin*latcos*rft; // shouldn't this be minus?
    if ( m==1 && n==3 )
        return std::sqrt(21.)/4.*latcos*(5.*latsin*latsin-1.)*rft; // shouldn't this be minus?
    if ( m==2 && n==2 )
        return -std::sqrt(15./2.)/2.*latcos*latcos*rft;
    if ( m==2 && n==3 )
        return std::sqrt(105./2.)/2.*latcos*latcos*latsin*rft;
    if ( m==3 && n==3 )
        return std::sqrt(35.)/4.*latcos*latcos*latcos*rft; // shouldn't this be minus?
    return 0;
}

void spectral_transform_grid_analytic(
        const size_t knsmax,  // truncation + 1 (in)
        const size_t knsmaxFT,// truncation + 1 for Fourier transformation (in)
        double n,             // total wave number (implemented so far for n<4
        double m,             // zonal wave number (implemented so far for m<4, m<n
        bool real,            // true: test real part, false: test imaginary part
        Grid grid,            // call with something like Grid("O32")
        double rspecg[],      // spectral data, size (knsmax+1)*knsmax (out)
        double rgp[])         // resulting grid point data (out)
{
    std::ostream& out = Log::info(); // just for debugging
    int N = (knsmax+2)*(knsmax+1)/2;
    for( int jm=0; jm<2*N; jm++) rspecg[jm] = 0.;
    int k = 0;
    for( int jm=0; jm<=knsmax; jm++ )
        for( int jn=jm; jn<=knsmax; jn++ ) {
            if( jm==m && jn==n ) {
                if( real ) {
                    rspecg[2*k] = 1.;
                } else {
                    rspecg[2*k+1] = 1.;
                }
            }
            k++;
        }

    for( int jm=0; jm<grid.size(); jm++) rgp[jm] = 0.;

    if( grid::StructuredGrid(grid) ) {
        grid::StructuredGrid g(grid);
        int idx = 0;
        for( size_t j=0; j<g.ny(); ++j ) {
            double lat = g.y(j);

            for( size_t i=0; i<g.nx(j); ++i ) {
                double lon = g.x(i,j);

                // compute spherical harmonics:
                rgp[idx] = sphericalharmonics_analytic(n, m, real, lon, lat);
                idx++;
            }
        }
    } else {
        int idx = 0;
        for( PointXY p: grid.xy()) {
            double lon = p.x(), lat = p.y();
            // compute spherical harmonics:
            rgp[idx++] = sphericalharmonics_analytic(n, m, real, lon, lat);
        }
    }
}

double spectral_transform_test(
        double trc,           // truncation
        double n,             // total wave number (implemented so far for n<4
        double m,             // zonal wave number (implemented so far for m<4, m<n
        bool real,            // true: test real part, false: test imaginary part
        Grid g,               // call with something like Grid("O32")
        bool pointwise)       // use point function for unstructured mesh for testing purposes
{
    std::ostream& out = Log::info();
    int N = (trc+2)*(trc+1)/2;
    auto *rspecg       = new double[2*N];
    auto *rgp          = new double[g.size()];
    auto *rgp_analytic = new double[g.size()];

    // compute analytic solution (this also initializes rspecg and needs to be done before the actual transform):
    spectral_transform_grid_analytic(trc, trc, n, m, real, g, rspecg, rgp_analytic);
    // perform spectral transform:
    spectral_transform_grid(trc, trc, g, rspecg, rgp, pointwise);

    double rms = 0.;
    for( int jm=0; jm<g.size(); jm++ ) rms += std::pow(rgp[jm]-rgp_analytic[jm],2);

    delete [] rspecg;
    delete [] rgp;
    delete [] rgp_analytic;

    rms = std::sqrt(rms/g.size());
    out << "m=" << m << " n=" << n << " real:" << real << " structured:" << grid::StructuredGrid(g) << " error:" << rms << std::endl;

    return rms;
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
  int trc = 1280; // truncation + 1
  int N = (trc+2)*(trc+1)/2;
  atlas::array::ArrayT<double> zlfpol_(N);
  atlas::array::ArrayView<double,1> zlfpol = make_view<double,1>(zlfpol_);

  double lat = 50.;
  lat = std::acos(0.99312859918509488)*util::Constants::radiansToDegrees();
  compute_legendre(trc, lat, zlfpol);
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

CASE( "test_transgeneral_point" )
{
  std::ostream& out = Log::info();
  Log::info() << "test_transgeneral_point" << std::endl;
  double tolerance = 1.e-15;
  // test spectral transform up to wave number 3 by comparing
  // the result with the analytically computed spherical harmonics

  Grid g = grid::UnstructuredGrid( new std::vector<PointXY>{
    {20., 50.},
    {-20., 30.},
    {-89., 179.},
    {70., -101.}
  });

  int trc = 1279; // truncation

  int n = 3, m = 3;
  double rms = 0.;
  for( int m=0; m<=3; m++ )
      for( int n=m; n<=3; n++ ) {
          rms = spectral_transform_test(trc, n, m, true, g, true);
          EXPECT( rms < tolerance );
          rms = spectral_transform_test(trc, n, m, false, g, true);
          EXPECT( rms < tolerance );
      }

}

//-----------------------------------------------------------------------------

CASE( "test_transgeneral_unstructured" )
{
  std::ostream& out = Log::info();
  Log::info() << "test_transgeneral_unstructured" << std::endl;
  double tolerance = 1.e-15;
  // test spectral transform up to wave number 3 by comparing
  // the result with the analytically computed spherical harmonics

  Grid g = grid::UnstructuredGrid( new std::vector<PointXY>{
    {20., 50.},
    {-20., 30.},
    {-89., 179.},
    {70., -101.}
  });

  int trc = 1279; // truncation

  int n = 3, m = 3;
  double rms = 0.;
  for( int m=0; m<=3; m++ )
      for( int n=m; n<=3; n++ ) {
          rms = spectral_transform_test(trc, n, m, true, g, false);
          EXPECT( rms < tolerance );
          rms = spectral_transform_test(trc, n, m, false, g, false);
          EXPECT( rms < tolerance );
      }

}

//-----------------------------------------------------------------------------

CASE( "test_transgeneral_structured" )
{
  std::ostream& out = Log::info();
  Log::info() << "test_transgeneral_structured" << std::endl;
  double tolerance = 1.e-15;
  // test spectral transform up to wave number 3 by comparing
  // the result with the analytically computed spherical harmonics

  std::string grid_uid("O10");
  grid::StructuredGrid g (grid_uid);

  int trc = 1279; // truncation

  int n = 3, m = 3;
  double rms = 0.;
  for( int m=0; m<=3; m++ )
      for( int n=m; n<=3; n++ ) {
          rms = spectral_transform_test(trc, n, m, true, g, false);
          EXPECT( rms < tolerance );
          rms = spectral_transform_test(trc, n, m, false, g, false);
          EXPECT( rms < tolerance );
      }

}

//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas


int main(int argc, char **argv) {
    atlas::test::AtlasTransEnvironment env( argc, argv );
    int i = run_tests ( argc, argv, false );
    std::cout << atlas::Trace::report() << std::endl;
    return i;
}
