/*
 * (C) Copyright 2013 ECMWF.
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
#include "atlas/trans/local/LegendrePolynomials.h"
#include "atlas/trans/local/LegendreTransforms.h"
#include "atlas/trans/local/FourierTransforms.h"
#include "atlas/array/MakeView.h"
#include "atlas/runtime/Trace.h"
#include "transi/trans.h"

#include "tests/AtlasTestEnvironment.h"

#include <iomanip>
#include <chrono>

using namespace eckit;

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

void compute_legendre(
        const size_t trc,                  // truncation (in)
        const double& lat,                 // latitude in radians (in)
        array::ArrayView<double,1>& zlfpol)// values of associated Legendre functions, size (trc+1)*trc/2 (out)
{
    trans::compute_legendre_polynomials(trc,lat,zlfpol.data());
}

//-----------------------------------------------------------------------------

void legendre_transform(
        const size_t trc,                    // truncation (in)
        const size_t trcFT,                  // truncation for Fourier transformation (in)
        array::ArrayView<double,1>& rlegReal,// values of associated Legendre functions, size (trc+1)*trc/2 (out)
        array::ArrayView<double,1>& rlegImag,// values of associated Legendre functions, size (trc+1)*trc/2 (out)
        const array::ArrayView<double,1>& zlfpol,  // values of associated Legendre functions, size (trc+1)*trc/2 (in)
        const double rspecg[])               // spectral data, size (trc+1)*trc (in)
{
    trans::invtrans_legendre( trc, trcFT, trc, zlfpol.data(), 1, rspecg, rlegReal.data(), rlegImag.data() );
}

//-----------------------------------------------------------------------------

double fourier_transform(
        const size_t trcFT,
        array::ArrayView<double,1>& rlegReal,// values of associated Legendre functions, size (trc+1)*trc/2 (out)
        array::ArrayView<double,1>& rlegImag,// values of associated Legendre functions, size (trc+1)*trc/2 (out)
        const double lon) // radians
{
  double gp[1];
  trans::invtrans_fourier( trcFT, lon, 1, rlegReal.data(), rlegImag.data(), gp );
  return gp[0];
}

//-----------------------------------------------------------------------------
// Routine to compute the spectral transform by using a local Fourier transformation
// for a single point
//
// Author:
// Andreas Mueller *ECMWF*
//
double spectral_transform_point(
        const size_t trc,                  // truncation (in)
        const size_t trcFT,                // truncation for Fourier transformation (in)
        const double lon,                  // longitude in radians (in)
        const double lat,                  // latitude in radians (in)
        const double rspecg[])             // spectral data, size (trc+1)*trc (in)
{
    int N = (trc+2)*(trc+1)/2;
    ATLAS_TRACE();
    atlas::array::ArrayT<double> zlfpol_(N);
    atlas::array::ArrayView<double,1> zlfpol = make_view<double,1>(zlfpol_);

    atlas::array::ArrayT<double> rlegReal_(trcFT+1);
    atlas::array::ArrayView<double,1> rlegReal = make_view<double,1>(rlegReal_);

    atlas::array::ArrayT<double> rlegImag_(trcFT+1);
    atlas::array::ArrayView<double,1> rlegImag = make_view<double,1>(rlegImag_);

    // Legendre transform:
    compute_legendre(trc, lat, zlfpol);
    legendre_transform(trc, trcFT, rlegReal, rlegImag, zlfpol, rspecg);

    // Fourier transform:
    return fourier_transform(trcFT, rlegReal, rlegImag, lon);
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
        const size_t trc,     // truncation (in)
        const size_t trcFT,   // truncation for Fourier transformation (in)
        const Grid grid,            // call with something like Grid("O32")
        const double rspecg[],      // spectral data, size (trc+1)*trc (in)
        double rgp[],         // resulting grid point data (out)
        const bool pointwise)       // use point function for unstructured mesh for testing purposes
{
    std::ostream& out = Log::info(); // just for debugging
    int N = (trc+2)*(trc+1)/2;
    ATLAS_TRACE();
    atlas::array::ArrayT<double> zlfpol_(N);
    atlas::array::ArrayView<double,1> zlfpol = make_view<double,1>(zlfpol_);

    atlas::array::ArrayT<double> rlegReal_(trcFT+1);
    atlas::array::ArrayView<double,1> rlegReal = make_view<double,1>(rlegReal_);

    atlas::array::ArrayT<double> rlegImag_(trcFT+1);
    atlas::array::ArrayView<double,1> rlegImag = make_view<double,1>(rlegImag_);

    int idx = 0;

    if( grid::StructuredGrid(grid) ) {
        grid::StructuredGrid g(grid);
        for( size_t j=0; j<g.ny(); ++j ) {
            double lat = g.y(j) * util::Constants::degreesToRadians() ;

            // Legendre transform:
            compute_legendre(trc, lat, zlfpol);
            legendre_transform(trc, trcFT, rlegReal, rlegImag, zlfpol, rspecg);

            for( size_t i=0; i<g.nx(j); ++i ) {
                double lon = g.x(i,j) * util::Constants::degreesToRadians();
                // Fourier transform:
                rgp[idx++] = fourier_transform(trcFT, rlegReal, rlegImag, lon);
            }
        }
    } else {
        for( PointXY p: grid.xy()) {
            double lon = p.x() * util::Constants::degreesToRadians();
            double lat = p.y() * util::Constants::degreesToRadians();
            if( pointwise ) {
                // alternative for testing: use spectral_transform_point function:
                rgp[idx++] = spectral_transform_point(trc, trcFT, lon, lat, rspecg);
            } else {
                // Legendre transform:
                compute_legendre(trc, lat, zlfpol);
                legendre_transform(trc, trcFT, rlegReal, rlegImag, zlfpol, rspecg);

                // Fourier transform:
                rgp[idx++] = fourier_transform(trcFT, rlegReal, rlegImag, lon);
            }
        }
    }

    EXPECT( idx == grid.size() );
}

//-----------------------------------------------------------------------------
// Routine to compute the spherical harmonics analytically at one point
// (up to wave number 3)
//
// Author:
// Andreas Mueller *ECMWF*
//
double sphericalharmonics_analytic_point(
        const double n,             // total wave number (implemented so far for n<4
        const double m,             // zonal wave number (implemented so far for m<4, m<n
        const int imag,             // 0: test real part, 1: test imaginary part
        const double lon,           // longitude in radians
        const double lat,           // latitude in radians
        const int ivar_in,          // variable that is set to 1 for wave number n,m. 0: vorticity, 1: divergence, 2: scalar
        const int ivar_out)         // variable returned by this function. 0: u, 1: v, 2: scalar
{
    double latsin = std::sin(  lat), latcos = std::cos(  lat);
    double lonsin = std::sin(m*lon), loncos = std::cos(m*lon);
    double a = util::Earth::radiusInMeters();
    // Fourier part of the spherical harmonics:
    double rft = 1.;
    if( m>0 ) rft *= 2.; // the famous factor 2 that noone really understands
    if( imag==0 ) {
        rft *= loncos;
    } else {
        rft *= -lonsin;
    }
    // Legendre part of the spherical harmonics (following http://mathworld.wolfram.com/SphericalHarmonic.html
    // multiplied with -2*sqrt(pi) due to different normalization and different coordinates):
    // (can also be computed on http://www.wolframalpha.com with:
    // LegendreP[n, m, x]/Sqrt[1/2*Integrate[LegendreP[n, m, y]^2, {y, -1, 1}]])
    // n, m need to be replaced by hand with the correct values
    // (otherwise the command will be too long for the free version of wolframalpha)

    // scalar:
    if ( ivar_in==2 ) {
        if ( ivar_out==2 ) {
            if ( m==0 && n==0 )
                return rft;
            if ( m==0 && n==1 )
                return std::sqrt(3.)*latsin*rft;
            if ( m==0 && n==2 )
                return std::sqrt(5.)/2.*(3.*latsin*latsin-1.)*rft; // sign?
            if ( m==0 && n==3 )
                return std::sqrt(7.)/2.*(5.*latsin*latsin-3.)*latsin*rft; // sign?
            if ( m==1 && n==1 )
                return std::sqrt(3./2.)*latcos*rft; // sign?
            if ( m==1 && n==2 )
                return std::sqrt(15./2.)*latsin*latcos*rft; // sign?
            if ( m==1 && n==3 )
                return std::sqrt(21.)/4.*latcos*(5.*latsin*latsin-1.)*rft; // sign?
            if ( m==2 && n==2 )
                return std::sqrt(15./2.)/2.*latcos*latcos*rft;
            if ( m==2 && n==3 )
                return std::sqrt(105./2.)/2.*latcos*latcos*latsin*rft;
            if ( m==3 && n==3 )
                return std::sqrt(35.)/4.*latcos*latcos*latcos*rft; // sign?
            if ( m==4 && n==4 )
                return (3*std::sqrt(17.5)*std::pow(latcos,4))/8.*rft;
            if ( m==5 && n==5 )
                return (3*std::sqrt(77)*std::pow(latcos,5))/16.*rft;
            if ( m==6 && n==6 )
                return (std::sqrt(3003)*std::pow(latcos,6))/32.*rft;
            if ( m==7 && n==7 )
                return (3*std::sqrt(357.5)*std::pow(latcos,7))/32.*rft;
            if ( m==8 && n==8 )
                return (3*std::sqrt(6077.5)*std::pow(latcos,8))/128.*rft;
            if ( m==9 && n==9 )
                return (std::sqrt(230945)*std::pow(latcos,9))/256.*rft;
            if ( m==10 && n==10 )
                return (std::sqrt(969969)*std::pow(latcos,10))/512.*rft;
            if ( m==11 && n==11 )
                return (std::sqrt(1.0140585e6)*std::pow(latcos,11))/512.*rft;
            if ( m==12 && n==12 )
                return (5*std::sqrt(676039)*std::pow(latcos,12))/2048.*rft;
            if ( m==13 && n==13 )
                return (15*std::sqrt(78004.5)*std::pow(latcos,13))/2048.*rft;
            if ( m==14 && n==14 )
                return (15*std::sqrt(323161.5)*std::pow(latcos,14))/4096.*rft;
            if ( m==15 && n==15 )
                return (3*std::sqrt(33393355)*std::pow(latcos,15))/8192.*rft;
            if ( m==16 && n==16 )
                return (3*std::sqrt(5.509903575e8)*std::pow(latcos,16))/32768.*rft;
            if ( m==17 && n==17 )
                return (15*std::sqrt(90751353)*std::pow(latcos,17))/65536.*rft;
            if ( m==18 && n==18 )
                return (5*std::sqrt(3357800061)*std::pow(latcos,18))/131072.*rft;
            if ( m==19 && n==19 )
                return (15*std::sqrt(3.829070245e8)*std::pow(latcos,19))/131072.*rft;
            if ( m==20 && n==20 )
                return (3*std::sqrt(156991880045)*std::pow(latcos,20))/524288.*rft;
            if ( m==21 && n==21 )
                return (std::sqrt(1.4465680375575e12)*std::pow(latcos,21))/524288.*rft;
            if ( m==22 && n==22 )
                return (15*std::sqrt(2.63012370465e10)*std::pow(latcos,22))/1.048576e6*rft;
            if ( m==23 && n==23 )
                return (15*std::sqrt(107492012277)*std::pow(latcos,23))/2.097152e6*rft;
            if ( m==24 && n==24 )
                return (105*std::sqrt(35830670759)*std::pow(latcos,24))/8.388608e6*rft;
            if ( m==25 && n==25 )
                return (21*std::sqrt(9.136821043545e11)*std::pow(latcos,25))/8.388608e6*rft;
            if ( m==26 && n==26 )
                return (21*std::sqrt(3.7250116562145e12)*std::pow(latcos,26))/1.6777216e7*rft;
            if ( m==27 && n==27 )
                return (7*std::sqrt(136583760727865.)*std::pow(latcos,27))/3.3554432e7*rft;
            if ( m==28 && n==28 )
                return (std::sqrt(2.7248460265209068e16)*std::pow(latcos,28))/6.7108864e7*rft;
            if ( m==29 && n==29 )
                return (std::sqrt(110873045217057585.)*std::pow(latcos,29))/1.34217728e8*rft;
            if ( m==30 && n==30 )
                return (std::sqrt(450883717216034179.)*std::pow(latcos,30))/2.68435456e8*rft;
            if ( m==31 && n==31 )
                return (21*std::sqrt(1.0389025742304935e15)*std::pow(latcos,31))/2.68435456e8*rft;
            if ( m==32 && n==32 )
                return (21*std::sqrt(6.752866732498208e16)*std::pow(latcos,32))/2.147483648e9*rft;
            if ( m==33 && n==33 )
                return (7*std::sqrt(2467865842240254105.)*std::pow(latcos,33))/4.294967296e9*rft;
            if ( m==34 && n==34 )
                return (21*std::sqrt(1112959105324036165.)*std::pow(latcos,34))/8.589934592e9*rft;
            if ( m==35 && n==35 )
                return (3*std::sqrt(5.53140675346046e19)*std::pow(latcos,35))/8.589934592e9*rft;
            if ( m==36 && n==36 )
                return (std::sqrt(8075853860052271220473.)*std::pow(latcos,36))/3.4359738368e10*rft;
            if ( m==37 && n==37 )
                return (5*std::sqrt(3.2739948081292994e20)*std::pow(latcos,37))/3.4359738368e10*rft;
            if ( m==38 && n==38 )
                return (35*std::sqrt(2.707815254843781e19)*std::pow(latcos,38))/6.8719476736e10*rft;
            if ( m==39 && n==39 )
                return (35*std::sqrt(109701233401363445369.)*std::pow(latcos,39))/1.37438953472e11*rft;
            if ( m==40 && n==40 )
                return (63*std::sqrt(548506167006817226845.)*std::pow(latcos,40))/5.49755813888e11*rft;
            if ( m==41 && n==41 )
                return (63*std::sqrt(5.551952666044613e20)*std::pow(latcos,41))/5.49755813888e11*rft;
            if ( m==42 && n==42 )
                return (15*std::sqrt(3.964094203555854e22)*std::pow(latcos,42))/1.099511627776e12*rft;
            if ( m==43 && n==43 )
                return (45*std::sqrt(17823059209786010066579.)*std::pow(latcos,43))/2.199023255552e12*rft;
            if ( m==44 && n==44 )
                return (45*std::sqrt(7.210237589413431e22)*std::pow(latcos,44))/4.398046511104e12*rft;
            if ( m==45 && n==45 )
                return (21*std::sqrt(1339044123748208678378695.)*std::pow(latcos,45))/8.796093022208e12*rft;
            //    return std::pow(latcos,45)*rft*21.*std::sqrt(1339044123748208678378695.)/8796093022208.; // sign?
        } else {
            return 0.;
        }
    }

    // for the remainder the factor 2 from rft is already included in the formulas:

    // vorticity:
    if ( ivar_in==0 ) {
        if ( ivar_out==0 ) { // u:
            if ( m==0 && n==0 )
                return 0.;
            if ( m==0 && n==1 ) {
                if ( imag==0 ) {
                    return std::sqrt(3.)*a/2.*latcos;
                } else {
                    return 0.;
                }
            }
            if ( m==1 && n==1 ) {
                if ( imag==0 ) {
                    return -a*std::sqrt(3./2.)*loncos*latsin;
                } else {
                    return a*std::sqrt(3./2.)*lonsin*latsin;
                }
            }
        } else if ( ivar_out==1 ) { // v:
            if ( m==0 && n==0 )
                return 0.;
            if ( m==0 && n==1 )
                return 0.;
            if ( m==1 && n==1 ) {
                if ( imag==0 ) {
                    return a*std::sqrt(3./2.)*lonsin;
                } else {
                    return a*std::sqrt(3./2.)*loncos;
                }
            }
        } else {
            return 0.;
        }
    }

    // divergence:
    if ( ivar_in==1 ) {
        if ( ivar_out==0 ) { // u:
            if ( m==0 && n==0 )
                return 0.;
            if ( m==0 && n==1 )
                return 0.;
            if ( m==1 && n==1 ) {
                if ( imag==0 ) {
                    return a*std::sqrt(3./2.)*lonsin;
                } else {
                    return a*std::sqrt(3./2.)*loncos;
                }
            }
        } else if ( ivar_out==1 ) { // v:
            if ( m==0 && n==0 )
                return 0.;
            if ( m==0 && n==1 ) {
                if ( imag==0 ) {
                    return -std::sqrt(3.)*a/2.*latcos;
                } else {
                    return 0.;
                }
            }
            if ( m==1 && n==1 ) {
                if ( imag==0 ) {
                    return a*std::sqrt(3./2.)*loncos*latsin;
                } else {
                    return -a*std::sqrt(3./2.)*lonsin*latsin; // sign?
                }
            }
        } else {
            return 0.;
        }
    }

    return -1.;
}

//-----------------------------------------------------------------------------
// Routine to compute the spherical harmonics analytically on a grid
// (up to wave number 3)
//
// Author:
// Andreas Mueller *ECMWF*
//
void spectral_transform_grid_analytic(
        const size_t trc,     // truncation (in)
        const size_t trcFT,   // truncation for Fourier transformation (in)
        const double n,       // total wave number (implemented so far for n<4
        const double m,       // zonal wave number (implemented so far for m<4, m<n
        const int imag,       // 0: test real part, 1: test imaginary part
        const Grid grid,      // call with something like Grid("O32")
        double rspecg[],      // spectral data, size (trc+1)*trc (out)
        double rgp[],         // resulting grid point data (out)
        const int ivar_in,    // variable that is set to 1 for wave number n,m. 0: vorticity, 1: divergence, 2: scalar
        const int ivar_out)   // variable returned by this function. 0: u, 1: v, 2: scalar
{
    int N = (trc+2)*(trc+1)/2;
    for( int jm=0; jm<2*N; jm++) rspecg[jm] = 0.;
    int k = 0;
    for( int jm=0; jm<=trc; jm++ ) {
        for( int jn=jm; jn<=trc; jn++ ) {
            if( jm==m && jn==n ) {
              rspecg[2*k+imag] = 1.;
              rspecg[2*k+(1-imag)] = 0.;
            }
            k++;
        }
    }

    for( int jm=0; jm<grid.size(); jm++) rgp[jm] = 0.;

    if( grid::StructuredGrid(grid) ) {
        grid::StructuredGrid g(grid);
        int idx = 0;
        bool fullgrid = (g.nx(0)==g.nxmax());
        for( size_t j=0; j<g.ny(); ++j ) {
            double lat = g.y(j) * util::Constants::degreesToRadians();

            for( size_t i=0; i<g.nx(j); ++i ) {
                double lon = g.x(i,j) * util::Constants::degreesToRadians();

                // compute spherical harmonics:
                if( trans::fourier_truncation(trc, g.nx(j), g.nxmax(), g.ny(), lat, fullgrid)>=m ) {
                    rgp[idx++] = sphericalharmonics_analytic_point(n, m, imag, lon, lat, ivar_in, ivar_out);
                } else {
                    rgp[idx++] = 0.;
                }
            }
        }
    } else {
        int idx = 0;
        for( PointXY p: grid.xy()) {
            double lon = p.x() * util::Constants::degreesToRadians();
            double lat = p.y() * util::Constants::degreesToRadians();
            // compute spherical harmonics:
            rgp[idx++] = sphericalharmonics_analytic_point(n, m, imag, lon, lat, ivar_in, ivar_out);
        }
    }
}

//-----------------------------------------------------------------------------
// Compute root mean square of difference between two arrays divided by maximum absolute value of second array
//
// Author:
// Andreas Mueller *ECMWF*
//
double compute_rms(
        const size_t N,  // length of the arrays
        double array1[], // first of the two arrays
        double array2[]) // second of the two arrays
{
    double rms = 0., rmax = 0.;
    for( int idx=0; idx<N; idx++ ) {
      double diff = array1[idx]-array2[idx];
      rms += diff*diff;
      rmax = std::max(rmax, std::abs(array2[idx]));
    }
    if( rmax==0. ) {
        rms = 0.;
    } else {
        rms = std::sqrt(rms/N)/rmax;
    }
    return rms;
}

//-----------------------------------------------------------------------------
// Routine to test the spectral transform by comparing it with the analytically
// derived spherical harmonics
//
// Author:
// Andreas Mueller *ECMWF*
//
double spectral_transform_test(
        double trc,           // truncation
        double n,             // total wave number (implemented so far for n<4
        double m,             // zonal wave number (implemented so far for m<4, m<n
        int imag,             // 0: test real part, 1: test imaginary part
        Grid g,               // call with something like Grid("O32")
        bool pointwise)       // use point function for unstructured mesh for testing purposes
{
    std::ostream& out = Log::info();
    int N = (trc+2)*(trc+1)/2;
    auto *rspecg       = new double[2*N];
    auto *rgp          = new double[g.size()];
    auto *rgp_analytic = new double[g.size()];

    // compute analytic solution (this also initializes rspecg and needs to be done before the actual transform):
    spectral_transform_grid_analytic(trc, trc, n, m, imag, g, rspecg, rgp_analytic, 2, 2);
    // perform spectral transform:

    spectral_transform_grid(trc, trc, g, rspecg, rgp, pointwise);

    //for( int i=0; i<g.size(); ++i ) rgp[i] = 0.;

    double rms = compute_rms(g.size(), rgp, rgp_analytic);

    /*out << "m=" << m << " n=" << n << " imag:" << imag << " structured:" << grid::StructuredGrid(g) << " error:" << rms;
    if( rms > 2.e-15 ) {
        out << " !!!!" << std::endl;
        for( int jp=0; jp<g.size(); jp++ ) {
            out << rgp[jp]/rgp_analytic[jp] << " rgp:" << rgp[jp] << " analytic:" << rgp_analytic[jp] << std::endl;
        }
    }
    out << std::endl;*/

    delete [] rspecg;
    delete [] rgp;
    delete [] rgp_analytic;

    return rms;
}

//-----------------------------------------------------------------------------
#if 0
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

  double lat = std::acos(0.99312859918509488);
  compute_legendre(trc, lat, zlfpol);
}
#endif
//-----------------------------------------------------------------------------
#if 1
CASE( "test_transgeneral_point" )
{
  std::ostream& out = Log::info();
  Log::info() << "test_transgeneral_point" << std::endl;
  double tolerance = 2.e-15;
  // test spectral transform up to wave number 3 by comparing
  // the result with the analytically computed spherical harmonics

  Grid g = grid::UnstructuredGrid( {
                                     {  50.,  20.},
                                     {  30., -20.},
                                     { 179., -89.},
                                     {-101.,  70.}
                                   } );

  int trc = 47; // truncation

  double rms = 0.;
  for( int m=0; m<=3; m++ ) { // zonal wavenumber
      for( int n=m; n<=3; n++ ) { // total wavenumber
          for( int imag=0; imag<=1; imag++ ) { // real and imaginary part
              rms = spectral_transform_test(trc, n, m, imag, g, true);
              EXPECT( rms < tolerance );
          }
      }
  }

}
#endif
//-----------------------------------------------------------------------------
#if 1
CASE( "test_transgeneral_unstructured" )
{
  std::ostream& out = Log::info();
  Log::info() << "test_transgeneral_unstructured" << std::endl;
  double tolerance = 2.e-15;
  // test spectral transform up to wave number 3 by comparing
  // the result with the analytically computed spherical harmonics

  Grid g = grid::UnstructuredGrid( new std::vector<PointXY>{
                                       {  50.,  20.},
                                       {  30., -20.},
                                       { 179., -89.},
                                       {-101.,  70.}
  });

  int trc = 47; // truncation

  double rms = 0.;
  for( int m=0; m<=3; m++ ) { // zonal wavenumber
      for( int n=m; n<=3; n++ ) { // total wavenumber
          for( int imag=0; imag<=1; imag++ ) { // real and imaginary part
              rms = spectral_transform_test(trc, n, m, imag, g, false);
              EXPECT( rms < tolerance );
          }
      }
  }

}

//-----------------------------------------------------------------------------

CASE( "test_transgeneral_structured" )
{
  std::ostream& out = Log::info();
  Log::info() << "test_transgeneral_structured" << std::endl;
  double tolerance = 2.e-15;
  // test spectral transform up to wave number 3 by comparing
  // the result with the analytically computed spherical harmonics

  std::string grid_uid("O10");
  grid::StructuredGrid g (grid_uid);

  int trc = 47; // truncation

  double rms = 0.;
  for( int m=0; m<=3; m++ ) { // zonal wavenumber
      for( int n=m; n<=3; n++ ) { // total wavenumber
          for( int imag=0; imag<=1; imag++ ) { // real and imaginary part
              rms = spectral_transform_test(trc, n, m, imag, g, false);
              EXPECT( rms < tolerance );
          }
      }
  }

}

//-----------------------------------------------------------------------------

CASE( "test_transgeneral_with_translib" )
{
  Log::info() << "test_transgeneral_with_translib" << std::endl;
  // test transgeneral by comparing its result with the trans library
  // this test is based on the test_nomesh case in test_trans.cc

  std::ostream& out = Log::info();
  double tolerance = 1.e-13;
  Grid g( "F24" );
  grid::StructuredGrid gs(g);
  int trc = 47;
  trans::Trans trans(g,trc) ;

  functionspace::Spectral          spectral   (trans);
  functionspace::StructuredColumns gridpoints (g);

  Field spf  = spectral.  createField<double>(option::name("spf"));
  Field gpf  = gridpoints.createField<double>(option::name("gpf"));

  int N = (trc+2)*(trc+1)/2;
  std::vector<double> rspecg       (2*N);
  std::vector<double> rgp          (g.size());
  std::vector<double> rgp_analytic (g.size());

  int k = 0;
  for( int m=0; m<=trc; m++ ) { // zonal wavenumber
      for( int n=m; n<=trc; n++ ) { // total wavenumber
          for( int imag=0; imag<=1; imag++ ) { // real and imaginary part

              if( sphericalharmonics_analytic_point(n, m, true, 0., 0., 2, 2) == 0. ) {

                  array::ArrayView<double,1> sp = array::make_view<double,1>(spf);
                  sp.assign(0.);
                  sp(k) = 1.;

                  EXPECT_NO_THROW( trans.invtrans(spf,gpf) );

                  spectral_transform_grid_analytic(trc, trc, n, m, imag, g, rspecg.data(), rgp_analytic.data(), 2, 2);

                  // compute spectral transform with the general transform:
                  spectral_transform_grid(trc, trc, g, sp.data(), rgp.data(), false);

                  array::ArrayView<double,1> gp = array::make_view<double,1>(gpf);
                  double rms_trans = compute_rms(g.size(), gp.data(), rgp.data());
                  double rms_gen   = compute_rms(g.size(), rgp.data(), rgp_analytic.data());

                  if( rms_gen >= tolerance ) {
                    ATLAS_DEBUG_VAR(rms_gen);
                    ATLAS_DEBUG_VAR(tolerance);
                  }
                  EXPECT( rms_gen < tolerance );
                  EXPECT( rms_trans < tolerance );
              }
              k++;
          }
      }
  }
}

//-----------------------------------------------------------------------------

CASE( "test_trans_vordiv_with_translib" )
{
  Log::info() << "test_trans_vordiv_with_translib" << std::endl;
  // test transgeneral by comparing its result with the trans library
  // this test is based on the test_nomesh case in test_trans.cc

  std::ostream& out = Log::info();
  double tolerance = 1.e-13;

  // resolution: (Reduce this number if the test takes too long!)
  int res = 12;
  
  Grid g( "O" + std::to_string(res) );
  grid::StructuredGrid gs(g);
  int trc = res*2-1;
  trans::Trans trans     (g, trc) ;
  trans::Trans transLocal(g, trc, util::Config("type","local"));

  functionspace::Spectral          spectral   (trans);
  functionspace::StructuredColumns gridpoints (g);

  int nb_scalar = 2, nb_vordiv = 2;
  int N = (trc+2)*(trc+1)/2, nb_all = nb_scalar+2*nb_vordiv;
  std::vector<double> sp           (2*N*nb_scalar);
  std::vector<double> vor          (2*N*nb_vordiv);
  std::vector<double> div          (2*N*nb_vordiv);
  std::vector<double> rspecg       (2*N);
  std::vector<double> gp           (nb_all*g.size());
  std::vector<double> rgp          (nb_all*g.size());
  std::vector<double> rgp_analytic (g.size());

  int icase = 0;
  for( int ivar_in=0; ivar_in<3; ivar_in++ ) { // vorticity, divergence, scalar
      for( int ivar_out=0; ivar_out<3; ivar_out++ ) { // u, v, scalar
          int nb_fld = 1;
          if( ivar_out==2) {
              tolerance = 1.e-13;
              nb_fld = nb_scalar;
          } else {
              tolerance = 2.e-6;
              nb_fld = nb_vordiv;
          }
          for( int jfld=0; jfld<nb_fld; jfld++ ) { // multiple fields
              int k = 0;
              for( int m=0; m<=trc; m++ ) { // zonal wavenumber
                  for( int n=m; n<=trc; n++ ) { // total wavenumber
                      for( int imag=0; imag<=1; imag++ ) { // real and imaginary part

                          if( sphericalharmonics_analytic_point(n, m, true, 0., 0., ivar_in, ivar_in) == 0. ) {

                              for( int j=0; j<2*N*nb_scalar; j++ ) {
                                  sp [j] = 0.;
                              }
                              for( int j=0; j<2*N*nb_vordiv; j++ ) {
                                  vor[j] = 0.;
                                  div[j] = 0.;
                              }
                              if( ivar_in==0 ) vor[k*nb_vordiv+jfld] = 1.;
                              if( ivar_in==1 ) div[k*nb_vordiv+jfld] = 1.;
                              if( ivar_in==2 ) sp [k*nb_scalar+jfld] = 1.;

                              for( int j=0; j<nb_all*g.size(); j++ ) {
                                  gp [j] = 0.;
                                  rgp[j] = 0.;
                              }
                              for( int j=0; j<g.size(); j++ ) {
                                  rgp_analytic[j] = 0.;
                              }

                              spectral_transform_grid_analytic(trc, trc, n, m, imag, g, rspecg.data(), rgp_analytic.data(), ivar_in, ivar_out);

                              EXPECT_NO_THROW( trans.invtrans( nb_scalar, sp.data(), nb_vordiv, vor.data(), div.data(), gp.data() ) );

                              // compute spectral transform with the general transform:
                              //EXPECT_NO_THROW( spectral_transform_grid(trc, trc, g, sp, rgp, false) );
                              //EXPECT_NO_THROW( transLocal.invtrans( nb_scalar, sp, rgp) );
                              EXPECT_NO_THROW( transLocal.invtrans( nb_scalar, sp.data(), nb_vordiv, vor.data(), div.data(), rgp.data()) );

                              int pos = (ivar_out*nb_vordiv+jfld);
                              //Log::info() << "Case " << icase << " Analytic solution: ivar_in=" << ivar_in << " ivar_out=" << ivar_out << " m=" << m << " n=" << n << " imag=" << imag << " k=" << k << std::endl;
                              //for( int j=0; j<g.size(); j++ ) Log::info() << std::setprecision(2) << rgp_analytic[j] << " ";
                              //Log::info() << std::endl;
                              //Log::info() << "Trans library: ivar_in=" << ivar_in << " ivar_out=" << ivar_out << " m=" << m << " n=" << n << " imag=" << imag << " k=" << k << std::endl;
                              //for( int j=pos*g.size(); j<(pos+1)*g.size(); j++ ) Log::info() << gp[j] << " ";
                              //Log::info() << std::endl;
                              //Log::info() << "Local transform: ivar_in=" << ivar_in << " ivar_out=" << ivar_out << " m=" << m << " n=" << n << " imag=" << imag << " k=" << k << std::endl;
                              //Log::info() << "pos=" << pos << " pos*g.size()=" << pos*g.size() << " (pos+1)*g.size()=" << (pos+1)*g.size() << std::endl;
                              //for( int j=pos*g.size(); j<(pos+1)*g.size(); j++ ) Log::info() << rgp[j] << " ";
                              //Log::info() << std::endl;
                              //Log::info() << std::endl;

                              double rms_trans = compute_rms(g.size(),  gp.data()+pos*g.size(), rgp_analytic.data());
                              double rms_gen   = compute_rms(g.size(), rgp.data()+pos*g.size(), rgp_analytic.data());
                              double rms_diff  = compute_rms(g.size(), rgp.data()+pos*g.size(), gp.data()+pos*g.size());

                              if( rms_gen>=tolerance || rms_trans>=tolerance || rms_diff>=tolerance ) {
                                Log::info() << "Case " << icase << " ivar_in=" << ivar_in << " ivar_out=" << ivar_out << " m=" << m << " n=" << n << " imag=" << imag << " k=" << k << std::endl;
                                ATLAS_DEBUG_VAR(rms_gen);
                                ATLAS_DEBUG_VAR(rms_trans);
                                ATLAS_DEBUG_VAR(rms_diff);
                                ATLAS_DEBUG_VAR(tolerance);
                              }
                              EXPECT( rms_trans < tolerance );
                              EXPECT( rms_gen < tolerance );
                              icase++;
                          }
                          k++;
                      }
                  }
              }
          }
      }
  }
  Log::info() << "Vordiv+scalar comparison with trans: all " << icase << " cases successfully passed!" << std::endl;
}

//-----------------------------------------------------------------------------

CASE( "test_trans_invtrans" ) {

  trans::Trans trans( Grid("O64"), 63, util::Config("type","local") );

  std::vector<double> rspec(trans.spectralCoefficients());
  std::vector<double> rgp(trans.grid().size());

  // TODO: rspec needs proper initial data

  trans.invtrans(1,rspec.data(),rgp.data());

}
#endif

//-----------------------------------------------------------------------------

#if 0
CASE( "test_trans_fourier_truncation" )
{
  Log::info() << "test_trans_fourier_truncation" << std::endl;
  // test transgeneral by comparing its result with the trans library
  // this test is based on the test_nomesh case in test_trans.cc

  Grid g( "F640" );
  grid::StructuredGrid gs(g);
  int ndgl = gs.ny();
  //int trc = 2*ndgl; // extreme high truncation (below linear)
  int trc = ndgl-1; // linear
  //int trc = 2./3.*ndgl-1; // quadratic
  //int trc = ndgl/2. -1; // cubic
  trans::Trans trans(g, trc) ;
  bool fullgrid = (gs.nx(0)==gs.nxmax());
  for( int j=0; j<gs.ny(); j+=80 ) {
      double lat = gs.y(j) * util::Constants::degreesToRadians();
      int trcFT = trans::fourier_truncation(trc, gs.nx(j), gs.nxmax(), gs.ny(), lat, fullgrid);
      //Log::info() << trcFT << "         " << gs.nx(j) << std::endl;
  }
  // TODO: create some real criterion to test fourier_truncation. So far only comparison with trans library through print statements.
}
#endif
//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas


int main(int argc, char **argv) {
  return atlas::test::run< atlas::test::AtlasTransEnvironment >( argc, argv );
}
