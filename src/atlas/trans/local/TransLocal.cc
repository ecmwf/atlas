/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "atlas/trans/local/TransLocal.h"

#include "atlas/array.h"
#include "atlas/runtime/ErrorHandling.h"
#include "atlas/runtime/Log.h"
#include "atlas/parallel/mpi/mpi.h"

namespace atlas {
namespace trans {

namespace {
static TransBuilderGrid<TransLocal> builder("local");
}

// --------------------------------------------------------------------------------------------------------------------
// Helper functions
// --------------------------------------------------------------------------------------------------------------------
namespace { // anonymous

size_t legendre_size( const size_t truncation ) {
  return (truncation+2)*(truncation+1)/2;
}

int fourier_truncation( double truncation, double lat ) {
  return truncation;
}

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
        const size_t trc,  // truncation (in)
        const double lat,  // latitude in radians (in)
        double zlfpol[] )  // values of associated Legendre functions, size (trc+1)*trc/2 (out)
{
    array::ArrayT<int> idxmn_(trc+1,trc+1);
    array::ArrayView<int,2> idxmn = array::make_view<int,2>(idxmn_);
    int j = 0;
    for( int jm=0; jm<=trc; ++jm ) {
        for( int jn=jm; jn<=trc; ++jn ) {
            idxmn(jm,jn) = j++;
        }
    }

    array::ArrayT<double> zfn_(trc+1,trc+1);
    array::ArrayView<double,2> zfn = array::make_view<double,2>(zfn_);

    int iodd;

    // Belousov, Swarztrauber use zfn(0,0)=std::sqrt(2.)
    // IFS normalisation chosen to be 0.5*Integral(Pnm**2) = 1
    zfn(0,0) = 2.;
    for( int jn=1; jn<=trc; ++jn ) {
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

    // --------------------
    // 1. First two columns
    // --------------------
    double zdlx1 = (M_PI_2-lat); // theta
    double zdlx = std::cos(zdlx1); // cos(theta)
    double zdlsita = std::sqrt(1.-zdlx*zdlx); // sin(theta) (this is how trans library does it)

    zlfpol[0] = 1.;
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
    for( int jn=2; jn<=trc; jn+=2 ) {
        double zdlk = 0.5*zfn(jn,0);
        double zdlldn = 0.0;
        // represented by only even k
        for (int jk=2; jk<=jn; jk+=2 ) {
            // normalised ordinary Legendre polynomial == \overbar{P_n}^0
            zdlk = zdlk + zfn(jn,jk)*std::cos(jk*zdlx1);
            // normalised associated Legendre polynomial == \overbar{P_n}^1
            zdlldn = zdlldn + 1./std::sqrt(jn*(jn+1.))*zfn(jn,jk)*jk*std::sin(jk*zdlx1);
        }
        zlfpol[idxmn(0,jn)] = zdlk;
        zlfpol[idxmn(1,jn)] = zdlldn;
    }

    // odd N
    for( int jn=1; jn<=trc; jn+=2 ) {
        double zdlk = 0.5*zfn(jn,0);
        double zdlldn = 0.0;
        // represented by only even k
        for (int jk=1; jk<=jn; jk+=2 ) {
            // normalised ordinary Legendre polynomial == \overbar{P_n}^0
            zdlk = zdlk + zfn(jn,jk)*std::cos(jk*zdlx1);
            // normalised associated Legendre polynomial == \overbar{P_n}^1
            zdlldn = zdlldn + 1./std::sqrt(jn*(jn+1.))*zfn(jn,jk)*jk*std::sin(jk*zdlx1);
        }
        zlfpol[idxmn(0,jn)] = zdlk;
        zlfpol[idxmn(1,jn)] = zdlldn;
    }

    // --------------------------------------------------------------
    // 2. Diagonal (the terms 0,0 and 1,1 have already been computed)
    //    Belousov, equation (23)
    // --------------------------------------------------------------

    double zdls = zdl1sita*std::numeric_limits<double>::min();
    for( int jn=2; jn<=trc; ++jn ) {
        zlfpol[idxmn(jn,jn)] = zlfpol[idxmn(jn-1,jn-1)]*zdlsita*std::sqrt((2.*jn+1.)/(2.*jn));
        if( std::abs(zlfpol[idxmn(jn,jn)]) < zdls ) zlfpol[idxmn(jn,jn)] = 0.0;
    }

    // ---------------------------------------------
    // 3. General recurrence (Belousov, equation 17)
    // ---------------------------------------------

    for( int jn=3; jn<=trc; ++jn ) {
        for( int jm=2; jm<jn; ++jm ) {
            zlfpol[idxmn(jm,jn)] =
                    std::sqrt(((2.*jn+1.)*(jn+jm-1.)*(jn+jm-3.))/
                              ((2.*jn-3.)*(jn+jm   )*(jn+jm-2.))) * zlfpol[idxmn(jm-2,jn-2)]
                  - std::sqrt(((2.*jn+1.)*(jn+jm-1.)*(jn-jm+1.))/
                              ((2.*jn-1.)*(jn+jm   )*(jn+jm-2.))) * zlfpol[idxmn(jm-2,jn-1)] * zdlx
                  + std::sqrt(((2.*jn+1.)*(jn-jm   )           )/
                              ((2.*jn-1.)*(jn+jm   )           )) * zlfpol[idxmn(jm,jn-1  )] * zdlx;
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
        const size_t trc,       // truncation (in)
        const size_t trcFT,     // truncation for Fourier transformation (in)
        const double zlfpol[],  // values of associated Legendre functions, size (trc+1)*trc/2 (in)
        const double rspecg[],  // spectral data, size (trc+1)*trc (in)
        double rlegReal[],      // values of associated Legendre functions, size (trc+1)*trc/2 (out)
        double rlegImag[] )     // values of associated Legendre functions, size (trc+1)*trc/2 (out)
{
    // Legendre transformation:
    int k = 0;
    for( int jm=0; jm<=trcFT; ++jm ) {
        rlegReal[jm] = 0.;
        rlegImag[jm] = 0.;
        for( int jn=jm; jn<=trc; ++jn, ++k ) {
            // not completely sure where this factor 2 comes from. One possible explanation:
            // normalization of trigonometric functions in the spherical harmonics
            // integral over square of trig function is 1 for m=0 and 0.5 (?) for m>0
            rlegReal[jm] += 2. * rspecg[2*k]   * zlfpol[k];
            rlegImag[jm] += 2. * rspecg[2*k+1] * zlfpol[k];
        }
    }
    // Undo factor 2 for (jm == 0)
    rlegReal[0] /= 2.;
    rlegImag[0] /= 2.;

}

//-----------------------------------------------------------------------------
// Routine to compute the local Fourier transformation
//
// Author:
// Andreas Mueller *ECMWF*
//
double fourier_transform(
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


} // namespace anonymous

// --------------------------------------------------------------------------------------------------------------------
// Class TransLocal
// --------------------------------------------------------------------------------------------------------------------

TransLocal::TransLocal( const Cache& cache, const Grid& grid, const long truncation, const eckit::Configuration& config ) :
    grid_(grid),
    truncation_(truncation),
    precompute_ ( config.getBool("precompute", true) )
{
    if ( precompute_ ) {
        if( grid::StructuredGrid g = grid_ ) {
            ATLAS_TRACE("Precompute legendre structured");
            size_t size(0);
            legendre_begin_.resize( g.ny() );
            for( size_t j=0; j<g.ny(); ++j ) {
                legendre_begin_[j] = size;
                size  += legendre_size( truncation_ );
            }
            legendre_.resize(size);

            for( size_t j=0; j<g.ny(); ++j ) {
              double lat = g.y(j) * util::Constants::degreesToRadians();
              compute_legendre( truncation_ , lat, legendre_data(j) );
            }
        } else {
            ATLAS_TRACE("Precompute legendre unstructured");
            size_t size(0);
            legendre_begin_.resize( grid_.size() );
            for( size_t j=0; j<grid_.size(); ++j ) {
                legendre_begin_[j] = size;
                size  += legendre_size( truncation_ );
            }
            legendre_.resize(size);
            int j(0);
            for( PointXY p: grid_.xy() ) {
                double lat = p.y() * util::Constants::degreesToRadians();
                compute_legendre( truncation_ , lat, legendre_data(j++) );
            }
        }
    }
}

// --------------------------------------------------------------------------------------------------------------------

TransLocal::TransLocal( const Grid& grid, const long truncation, const eckit::Configuration& config ) :
  TransLocal( Cache(), grid, truncation, config ) {
}

// --------------------------------------------------------------------------------------------------------------------

TransLocal::~TransLocal()
{
}

// --------------------------------------------------------------------------------------------------------------------


void TransLocal::invtrans(
    const Field& spfield,
          Field& gpfield,
    const eckit::Configuration& config ) const
{
  NOTIMP;
}

// --------------------------------------------------------------------------------------------------------------------

void TransLocal::invtrans(
    const FieldSet& spfields,
          FieldSet& gpfields,
    const eckit::Configuration& config ) const
{
  NOTIMP;
}

// --------------------------------------------------------------------------------------------------------------------

void TransLocal::invtrans_grad(
    const Field& spfield,
          Field& gradfield,
    const eckit::Configuration& config ) const
{
  NOTIMP;
}

// --------------------------------------------------------------------------------------------------------------------

void TransLocal::invtrans_grad(
    const FieldSet& spfields,
          FieldSet& gradfields,
    const eckit::Configuration& config ) const
{
  NOTIMP;
}

// --------------------------------------------------------------------------------------------------------------------

void TransLocal::invtrans_vordiv2wind(
    const Field& spvor, const Field& spdiv,
    Field& gpwind,
    const eckit::Configuration& config ) const
{
  NOTIMP;
}

//-----------------------------------------------------------------------------
// Routine to compute the spectral transform by using a local Fourier transformation
// for a grid (same latitude for all longitudes, allows to compute Legendre functions
// once for all longitudes)
//
// Author:
// Andreas Mueller *ECMWF*
//
void TransLocal::invtrans(
    const int nb_scalar_fields, const double scalar_spectra[],
    double gp_fields[],
    const eckit::Configuration& config ) const
{
    ASSERT( nb_scalar_fields == 1 ); // More not yet supported

    // Depending on "precompute_legendre_", we have to compute the
    // legendre polynomials for every latitute
    std::vector<double> recomputed_legendre_;

    auto legPol = [&](double lat, int j) {
      if( precompute_ ) {
        return legendre_data(j);
      } else {
        recomputed_legendre_.resize( legendre_size( truncation_ ) );
        compute_legendre( truncation_, lat, recomputed_legendre_.data() );
        return recomputed_legendre_.data();
      }
    };

    // Temporary storage for legendre space
    std::vector<double> legReal(truncation_+1);
    std::vector<double> legImag(truncation_+1);

    // Transform
    if( grid::StructuredGrid g = grid_ ) {
        ATLAS_TRACE( "invtrans structured");
        int idx = 0;
        for( size_t j=0; j<g.ny(); ++j ) {
            double lat = g.y(j) * util::Constants::degreesToRadians();
            double trcFT = fourier_truncation( truncation_, lat );

            // Legendre transform:
            legendre_transform( truncation_, trcFT, legPol(lat,j), scalar_spectra, legReal.data(), legImag.data() );

            // Fourier transform:
            for( size_t i=0; i<g.nx(j); ++i ) {
                double lon = g.x(i,j) * util::Constants::degreesToRadians();
                gp_fields[idx++] = fourier_transform( trcFT, legReal.data(), legImag.data(), lon );
            }
        }
    } else {
        ATLAS_TRACE( "invtrans unstructured");
        int idx = 0;
        for( PointXY p: grid_.xy() ) {
            double lon = p.x() * util::Constants::degreesToRadians();
            double lat = p.y() * util::Constants::degreesToRadians();
            double trcFT = fourier_truncation( truncation_, lat );

            // Legendre transform:
            legendre_transform( truncation_, trcFT, legPol(lat,idx), scalar_spectra, legReal.data(), legImag.data() );

            // Fourier transform:
            gp_fields[idx++] = fourier_transform( trcFT, legReal.data(), legImag.data(), lon );
        }
    }
}

// --------------------------------------------------------------------------------------------------------------------

void TransLocal::invtrans(
    const int nb_vordiv_fields,
    const double vorticity_spectra[],
    const double divergence_spectra[],
    double gp_fields[],
    const eckit::Configuration& config ) const
{
  NOTIMP;
}

// --------------------------------------------------------------------------------------------------------------------

void TransLocal::invtrans(
    const int nb_scalar_fields,
    const double scalar_spectra[],
    const int nb_vordiv_fields,
    const double vorticity_spectra[],
    const double divergence_spectra[],
    double gp_fields[],
    const eckit::Configuration& config ) const
{
  NOTIMP;
}

// --------------------------------------------------------------------------------------------------------------------

void TransLocal::dirtrans(
    const Field& gpfield,
          Field& spfield,
    const eckit::Configuration& config ) const
{
  NOTIMP;
  // Not implemented and not planned.
  // Use the TransIFS implementation instead.
}

// --------------------------------------------------------------------------------------------------------------------

void TransLocal::dirtrans(
    const FieldSet& gpfields,
          FieldSet& spfields,
    const eckit::Configuration& config ) const
{
  NOTIMP;
  // Not implemented and not planned.
  // Use the TransIFS implementation instead.
}

// --------------------------------------------------------------------------------------------------------------------

void TransLocal::dirtrans_wind2vordiv(
    const Field& gpwind,
    Field& spvor, Field& spdiv,
    const eckit::Configuration& config ) const
{
  NOTIMP;
  // Not implemented and not planned.
  // Use the TransIFS implementation instead.
}

// --------------------------------------------------------------------------------------------------------------------

void TransLocal::dirtrans(
    const int nb_fields,
    const double scalar_fields[],
    double scalar_spectra[],
    const eckit::Configuration& ) const
{
  NOTIMP;
  // Not implemented and not planned.
  // Use the TransIFS implementation instead.
}

// --------------------------------------------------------------------------------------------------------------------

void TransLocal::dirtrans(
    const int nb_fields,
    const double wind_fields[],
    double vorticity_spectra[],
    double divergence_spectra[],
    const eckit::Configuration& ) const
{
  NOTIMP;
  // Not implemented and not planned.
  // Use the TransIFS implementation instead.
}

// --------------------------------------------------------------------------------------------------------------------

} // namespace trans
} // namespace atlas
