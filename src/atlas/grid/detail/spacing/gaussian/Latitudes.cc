/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

/// @file   Latitidus.cc
/// @author Willem Deconinck
/// @date   Jan 2014

#include <limits>
#include <memory>

#include "eckit/log/Bytes.h"

#include "atlas/array.h"
#include "atlas/grid/detail/spacing/gaussian/Latitudes.h"
#include "atlas/grid/detail/spacing/gaussian/N.h"
#include "atlas/library/config.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"
#include "atlas/runtime/Trace.h"
#include "atlas/util/Constants.h"
#include "atlas/util/CoordinateEnums.h"

//using eckit::ConcreteBuilderT0;
//using eckit::Factory;

using atlas::array::Array;
using atlas::array::ArrayView;
using atlas::array::make_view;

namespace atlas {
namespace grid {
namespace spacing {
namespace gaussian {

//-----------------------------------------------------------------------------

void compute_gaussian_quadrature_npole_equator( const size_t N, double lat[], double weights[] );

//-----------------------------------------------------------------------------

void gaussian_latitudes_npole_equator( const size_t N, double lats[] ) {
    std::stringstream Nstream;
    Nstream << N;
    std::string Nstr = Nstream.str();
    if ( GaussianLatitudesFactory::has( Nstr ) ) {
        std::unique_ptr<const GaussianLatitudes> gl( GaussianLatitudesFactory::build( Nstr ) );
        gl->assign( lats, N );
    }
    //    if ( Factory<GaussianLatitudes>::instance().exists( Nstr ) ) {
    //        std::unique_ptr<GaussianLatitudes> gl( Factory<GaussianLatitudes>::instance().get( Nstr ).create() );
    //        gl->assign( lats, N );
    //    }
    else {
        std::vector<double> weights( N );
        compute_gaussian_quadrature_npole_equator( N, lats, weights.data() );
    }
}

//-----------------------------------------------------------------------------

void gaussian_latitudes_npole_spole( const size_t N, double lats[] ) {
    gaussian_latitudes_npole_equator( N, lats );
    size_t end = 2 * N - 1;
    for ( size_t jlat = 0; jlat < N; ++jlat ) {
        lats[end--] = -lats[jlat];
    }
}

//-----------------------------------------------------------------------------

void gaussian_quadrature_npole_equator( const size_t N, double lats[], double weights[] ) {
    compute_gaussian_quadrature_npole_equator( N, lats, weights );
}

//-----------------------------------------------------------------------------

void gaussian_quadrature_npole_spole( const size_t N, double lats[], double weights[] ) {
    gaussian_quadrature_npole_equator( N, lats, weights );
    size_t end = 2 * N - 1;
    for ( size_t jlat = 0; jlat < N; ++jlat ) {
        lats[end]    = -lats[jlat];
        weights[end] = weights[jlat];
        end--;
    }
}

//-----------------------------------------------------------------------------

namespace {  // Anonymous namespace

//-----------------------------------------------------------------------------

void legpol_newton_iteration( int kn, const double pfn[], double px, double& pxn, double& pxmod ) {
    // Routine to perform a single Newton iteration step to find
    //                the zero of the ordinary Legendre polynomial of degree N

    //        Explicit arguments :
    //        --------------------
    //          KN       :  Degree of the Legendre polynomial              (in)
    //          KODD     :  odd or even number of latitudes                (in)
    //          PFN      :  Fourier coefficients of series expansion       (in)
    //                      for the ordinary Legendre polynomials
    //          PX       :  abcissa where the computations are performed   (in)
    //          PXN      :  new abscissa (Newton iteration)                (out)
    //          PXMOD    :  PXN-PX                                         (out)

    double zdlx, zdlk, zdlldn, zdlxn, zdlmod;
    int ik;
    int kodd = kn % 2;  // mod(kn,2)

    zdlx = px;
    zdlk = 0.;
    if ( kodd == 0 ) {
        zdlk = 0.5 * pfn[0];
    }
    zdlxn  = 0.;
    zdlldn = 0.;
    ik     = 1;

    // Newton interation step
    // ----------------------
    for ( int jn = 2 - kodd; jn <= kn; jn += 2 ) {
        // normalised ordinary Legendre polynomial == \overbar{P_n}^0
        zdlk += pfn[ik] * std::cos( static_cast<double>( jn ) * zdlx );
        // normalised derivative == d/d\theta(\overbar{P_n}^0)
        zdlldn -= pfn[ik] * static_cast<double>( jn ) * std::sin( static_cast<double>( jn ) * zdlx );
        ++ik;
    }
    // Newton method
    zdlmod = -zdlk / zdlldn;
    zdlxn  = zdlx + zdlmod;
    pxn    = zdlxn;
    pxmod  = zdlmod;
}

void legpol_weight( const int kn, const double pfn[], const double px, double& pw ) {
    // Routine to compute the quadrature weight of the legendre polynomial

    //        Explicit arguments :
    //        --------------------
    //          KN       :  Degree of the Legendre polynomial              (in)
    //          KODD     :  odd or even number of latitudes                (in)
    //          PFN      :  Fourier coefficients of series expansion       (in)
    //                      for the ordinary Legendre polynomials
    //          PX       :  abcissa where the computations are performed   (in)
    //          KFLAG    :  When KFLAG.EQ.1 computes the weights           (in)
    //          PXN      :  new abscissa (Newton iteration)                (out)
    //          PXMOD    :  PXN-PX                                         (out)

    double zdlx, zdlldn;
    int ik;
    int kodd = kn % 2;

    zdlx   = px;
    zdlldn = 0.;
    ik     = 1;

    // Compute weights
    // ---------------
    for ( int jn = 2 - kodd; jn <= kn; jn += 2 ) {
        zdlldn -= pfn[ik] * static_cast<double>( jn ) * std::sin( static_cast<double>( jn ) * zdlx );
        ++ik;
    }
    pw = static_cast<double>( 2 * kn + 1 ) / ( zdlldn * zdlldn );
}

//-----------------------------------------------------------------------------

void legpol_quadrature( const int kn, const double pfn[], double& pl, double& pw, int& kiter, double& pmod ) {
    //**** *GAWL * - Routine to perform the Newton loop

    //     Purpose.
    //     --------
    //           Find 0 of Legendre polynomial with Newton loop
    //**   Interface.
    //     ----------
    //        *CALL* *GAWL(PFN,PL,PW,PEPS,KN,KITER,PMOD)

    //        Explicit arguments :
    //        --------------------
    // PFN    Fourier coefficients of series expansion
    //        for the ordinary Legendre polynomials     (in)
    // PL     Gaussian latitude                         (inout)
    // PEPS   0 of the machine                          (in)
    // KN     Truncation                                (in)
    // KITER  Number of iterations                      (out)
    // PMOD   Last modification                         (inout)

    int iflag, itemax;

    double zx  = 0;
    double zw  = 0;
    double zxn = 0;

    //*       1. Initialization.
    //           ---------------

    itemax = 20;
    zx     = pl;
    iflag  = 0;
    pw     = 0.;

    //*       2. Newton iteration.
    //           -----------------

    const double zeps = std::numeric_limits<double>::epsilon();

    for ( int jter = 1; jter <= itemax + 1; ++jter ) {
        kiter = jter;
        legpol_newton_iteration( kn, pfn, zx, zxn, pmod );
        zx = zxn;
        if ( iflag == 1 ) {
            legpol_weight( kn, pfn, zx, zw );
            break;
        }
        if ( std::abs( pmod ) <= zeps * 1000. ) {
            iflag = 1;
        }
    }
    if ( iflag != 1 ) {
        std::stringstream s;
        s << "Could not converge gaussian latitude to accuracy [" << zeps * 1000 << "]\n";
        s << "after " << itemax << " iterations. Consequently also failed to compute quadrature weight.";
        throw_Exception( s.str(), Here() );
    }

    pl = zxn;
    pw = zw;
}

//-----------------------------------------------------------------------------

}  //  anonymous namespace

//-----------------------------------------------------------------------------

void compute_gaussian_quadrature_npole_equator( const size_t N, double lats[], double weights[] ) {
    Log::debug() << "Atlas computing Gaussian latitudes for N " << N << " which requires temporary memory of "
                 << eckit::Bytes( sizeof( double ) * ( 2 * N + 1 ) * ( 2 * N + 1 ) ) << std::endl;
    ATLAS_TRACE();

    // Compute first guess for colatitudes in radians
    double z;
    for ( size_t i = 0; i < N; ++i ) {
        z       = ( 4. * ( i + 1. ) - 1. ) * M_PI / ( 4. * 2. * N + 2. );
        lats[i] = ( z + 1. / ( tan( z ) * ( 8. * ( 2. * N ) * ( 2. * N ) ) ) );
    }

    int kdgl = 2 * N;
    array::ArrayT<double> zfn_( kdgl + 1, kdgl + 1 );
    //    WARNING: potentially HUGE allocation ( N=16000 --> 7.6 Gbytes )

    ArrayView<double, 2> zfn = make_view<double, 2>( zfn_ );

    int iodd;

    // Belousov, Swarztrauber use zfn(0,0)=std::sqrt(2.)
    // IFS normalisation chosen to be 0.5*Integral(Pnm**2) = 1
    zfn( 0, 0 ) = 2.;
    for ( int jn = 1; jn <= kdgl; ++jn ) {
        double zfnn = zfn( 0, 0 );
        for ( int jgl = 1; jgl <= jn; ++jgl ) {
            zfnn *= std::sqrt( 1. - 0.25 / ( static_cast<double>( jgl * jgl ) ) );
        }
        iodd          = jn % 2;
        zfn( jn, jn ) = zfnn;
        for ( int jgl = 2; jgl <= jn - iodd; jgl += 2 ) {
            zfn( jn, jn - jgl ) = zfn( jn, jn - jgl + 2 ) * static_cast<double>( ( jgl - 1 ) * ( 2 * jn - jgl + 2 ) ) /
                                  static_cast<double>( jgl * ( 2 * jn - jgl + 1 ) );
        }
    }

    iodd   = kdgl % 2;
    int ik = iodd;

    std::vector<double> zzfn( N + 1 );
    std::vector<double> zmod( kdgl );
    std::vector<int> iter( kdgl );

    for ( int jgl = iodd; jgl <= kdgl; jgl += 2 ) {
        zzfn[ik] = zfn( kdgl, jgl );
        ++ik;
    }

    const double pole = 90.;
    for ( size_t jgl = 0; jgl < N; ++jgl ) {
        // refine colat first guess here via Newton's method
        legpol_quadrature( kdgl, zzfn.data(), lats[jgl], weights[jgl], iter[jgl], zmod[jgl] );

        // Convert colat to lat, in degrees
        lats[jgl] = pole - lats[jgl] * util::Constants::radiansToDegrees();
    }
}

//-----------------------------------------------------------------------------

}  // namespace gaussian
}  // namespace spacing
}  // namespace grid
}  // namespace atlas
