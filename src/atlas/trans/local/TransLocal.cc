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
#include "atlas/trans/local/FourierTransforms.h"
#include "atlas/trans/local/LegendreTransforms.h"
#include "atlas/trans/local/LegendrePolynomials.h"
#include "atlas/util/Earth.h"

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
              compute_legendre_polynomials( truncation_ , lat, legendre_data(j) );
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
                compute_legendre_polynomials( truncation_ , lat, legendre_data(j++) );
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

void TransLocal::invtrans(
    const int nb_scalar_fields, const double scalar_spectra[],
    double gp_fields[],
    const eckit::Configuration& config ) const
{
    invtrans_uv(nb_scalar_fields, 0, scalar_spectra, gp_fields, config);
}

//-----------------------------------------------------------------------------
// Routine to compute the spectral transform by using a local Fourier transformation
// for a grid (same latitude for all longitudes, allows to compute Legendre functions
// once for all longitudes). U and v components are divided by cos(latitude) for
// nb_vordiv_fields > 0.
//
// Author:
// Andreas Mueller *ECMWF*
//
void TransLocal::invtrans_uv(
    const int nb_scalar_fields, const int nb_vordiv_fields, const double scalar_spectra[],
    double gp_fields[],
    const eckit::Configuration& config ) const
{
    int nb_fields = nb_scalar_fields;

    // Depending on "precompute_legendre_", we have to compute the
    // legendre polynomials for every latitute
    std::vector<double> recomputed_legendre_;

    auto legPol = [&](double lat, int j) -> const double * {
      if( precompute_ ) {
        return legendre_data(j);
      } else {
        recomputed_legendre_.resize( legendre_size( truncation_ ) );
        compute_legendre_polynomials( truncation_, lat, recomputed_legendre_.data() );
        return recomputed_legendre_.data();
      }
    };

    // Temporary storage for legendre space
    std::vector<double> legReal(nb_fields*(truncation_+1));
    std::vector<double> legImag(nb_fields*(truncation_+1));

    // Transform
    if( grid::StructuredGrid g = grid_ ) {
        ATLAS_TRACE( "invtrans structured");
        int idx = 0;
        for( size_t j=0; j<g.ny(); ++j ) {
            double lat = g.y(j) * util::Constants::degreesToRadians();
            double trcFT = fourier_truncation( truncation_, lat );

            // Legendre transform:
            invtrans_legendre( truncation_, trcFT, legPol(lat,j), nb_fields, scalar_spectra, legReal.data(), legImag.data() );

            // Fourier transform:
            for( size_t i=0; i<g.nx(j); ++i ) {
                double lon = g.x(i,j) * util::Constants::degreesToRadians();
                invtrans_fourier( trcFT, lon, nb_fields, legReal.data(), legImag.data(), gp_fields+(nb_fields*idx));
                ++idx;
            }

            // Divide U, V by cos(latitude):
            for( int jfld=2*nb_vordiv_fields; jfld<4*nb_vordiv_fields; ++jfld ) {
                gp_fields[nb_fields*idx+jfld] /= std::cos(lat);
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
            invtrans_legendre( truncation_, trcFT, legPol(lat,idx), nb_fields, scalar_spectra, legReal.data(), legImag.data() );

            // Fourier transform:
            invtrans_fourier( trcFT, lon, nb_fields, legReal.data(), legImag.data(), gp_fields+(nb_fields*idx));
            ++idx;

            // Divide U, V by cos(latitude):
            for( int jfld=2*nb_vordiv_fields; jfld<4*nb_vordiv_fields; ++jfld ) {
                gp_fields[nb_fields*idx+jfld] /= std::cos(lat);
            }
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
void prfi1b(
    const int truncation,
    const int km, // zonal wavenumber
    const int nb_fields, // number of fields
    const double rspec[], // spectral data
    double pia[]) // spectral components in data layout of trans library
{
    int ilcm = truncation-km, ioff = (2*truncation-km+1)*km, nlei1 = truncation+4+(truncation+4+1)%2;
    for( int j=0; j<=ilcm; ++j ) {
        int inm = ioff+(ilcm-j)*2;
        for( int jfld=0; jfld<nb_fields; ++jfld ) {
            int ir = 2*jfld+1, ii = ir+1;
            pia[ir*nlei1+j+2] = rspec[inm*nb_fields+jfld];
            pia[ii*nlei1+j+2] = rspec[inm*nb_fields+jfld+1];
        }
    }

    for( int jfld=0; jfld<2*nb_fields; ++jfld ) {
        pia[jfld*nlei1  ]      = 0.;
        pia[jfld*nlei1+1]      = 0.;
        pia[jfld*nlei1+ilcm+3] = 0.;
    }
}

// --------------------------------------------------------------------------------------------------------------------
void vd2uv(
    const int truncation,
    const int km, // zonal wavenumber
    const int nb_vordiv_fields,
    const double vorticity_spectra[],
    const double divergence_spectra[],
    int nb_all_fields,
    double all_spectra[],
    int nb_scalar_fields,
    const double scalar_spectra[],
    const eckit::Configuration& config )
{
    int nlei1 = truncation+4+(truncation+4+1)%2;

    double repsnm[(truncation+1)*(truncation+6)/2], rlapin[truncation+3];
    int idx = 0;
    for( int jm=0; jm<=truncation; ++jm ) {
        for( int jn=jm; jn<=truncation+2; ++jn, ++idx ) {
            repsnm[idx] = std::sqrt((jn*jn-jm*jm)/(4.*jn*jn-1.));
        }
    }
    for( int jn=0; jn<=truncation+1; ++jn ) {
        rlapin[jn] = -util::Earth::radiusInMeters()*util::Earth::radiusInMeters()/(jn*(jn+1.));
    }
    double rn[truncation+3];
    for( int jn=0; jn<truncation+3; ++jn ) {
        rn[jn] = jn;
    }
    double zepsnm[truncation+6], zlapin[truncation+6], zn[truncation+6];
    for( int jn=km-1; jn<=truncation+2; ++jn ) {
        int ij = truncation+3-jn;
        zlapin[ij] = rlapin[jn];
        zepsnm[ij] = repsnm[jn];
        zn[ij]     = jn;
    }
    zn[0] = truncation+3;

    double rvor[2*nb_vordiv_fields*nlei1], rdiv[2*nb_vordiv_fields*nlei1];
    double ru[2*nb_vordiv_fields*nlei1], rv[2*nb_vordiv_fields*nlei1];
    prfi1b(truncation, km, nb_vordiv_fields, vorticity_spectra,  rvor);
    prfi1b(truncation, km, nb_vordiv_fields, divergence_spectra, rdiv);

    if( km==0 ) {
        for( int jfld=0; jfld<=nb_vordiv_fields; ++jfld ) {
            int ir=2*jfld*nlei1;
            for( int ji=1; ji<truncation+3-km; ++ji ) {
                ru[ir+ji] = + zn[ji+1]*repsnm[ji  ]*rlapin[ji+1]*rvor[ir+ji+1]
                            - zn[ji-2]*repsnm[ji-1]*rlapin[ji-1]*rvor[ir+ji-1];
                rv[ir+ji] = - zn[ji+1]*repsnm[ji  ]*rlapin[ji+1]*rdiv[ir+ji+1]
                            + zn[ji-2]*repsnm[ji-1]*rlapin[ji-1]*rdiv[ir+ji-1];
            }
        }
    } else {
        for( int jfld=0; jfld<=nb_vordiv_fields; ++jfld) {
            int ir=2*jfld*nlei1, ii = ir+nlei1;
            for( int ji=1; ji<truncation+3-km; ++ji ) {
                ru[ir+ji] = -                    km*rlapin[ji  ]*rdiv[ii+ji  ]
                            + zn[ji+1]*repsnm[ji  ]*rlapin[ji+1]*rvor[ir+ji+1]
                            - zn[ji-2]*repsnm[ji-1]*rlapin[ji-1]*rvor[ir+ji-1];
                ru[ii+ji] = +                    km*rlapin[ji  ]*rdiv[ir+ji  ]
                            + zn[ji+1]*repsnm[ji  ]*rlapin[ji+1]*rvor[ii+ji+1]
                            - zn[ji-2]*repsnm[ji-1]*rlapin[ji-1]*rvor[ii+ji-1];
                rv[ir+ji] = -                    km*rlapin[ji  ]*rvor[ii+ji  ]
                            - zn[ji+1]*repsnm[ji  ]*rlapin[ji+1]*rdiv[ir+ji+1]
                            + zn[ji-2]*repsnm[ji-1]*rlapin[ji-1]*rdiv[ir+ji-1];
                rv[ii+ji] = -                    km*rlapin[ji  ]*rvor[ir+ji  ]
                            - zn[ji+1]*repsnm[ji  ]*rlapin[ji+1]*rdiv[ii+ji+1]
                            + zn[ji-2]*repsnm[ji-1]*rlapin[ji-1]*rdiv[ii+ji-1];
            }
        }
    }

    int ilcm = truncation-km, ioff = (2*truncation-km+1)*km;
    double za_r = 1./util::Earth::radiusInMeters();
    for( int j=0; j<ilcm; ++j ) {
        int inm = ioff+(ilcm-j)*2;
        for( int jfld=0; jfld<=nb_vordiv_fields; ++jfld ) { // < instead of <= ????
            int ir = 2*jfld*nlei1, ii = ir + nlei1;
            int idx = inm*nb_all_fields+jfld, idx0 = inm*nb_vordiv_fields+jfld, idx1 = idx0 + nb_vordiv_fields;
            // vorticity
            all_spectra[idx] = vorticity_spectra[idx0];
            idx += nb_vordiv_fields;
            all_spectra[idx] = vorticity_spectra[idx1];
            idx += nb_vordiv_fields;
            // divergence
            all_spectra[idx] = divergence_spectra[idx0];
            idx += nb_vordiv_fields;
            all_spectra[idx] = divergence_spectra[idx1];
            idx += nb_vordiv_fields;
            // u
            all_spectra[idx] = ru[ir+j+2]*za_r;
            idx += nb_vordiv_fields;
            all_spectra[idx] = ru[ii+j+2]*za_r;
            idx += nb_vordiv_fields;
            // v
            all_spectra[idx] = rv[ir+j+2]*za_r;
            idx += nb_vordiv_fields;
            all_spectra[idx] = rv[ii+j+2]*za_r;
        }
        for( int jfld=0; jfld<=nb_scalar_fields; ++jfld ) {
            int idx = inm*nb_all_fields+8*nb_vordiv_fields+jfld, idx0 = inm*nb_scalar_fields+jfld, idx1 = idx0 + nb_scalar_fields;
            // scalars
            all_spectra[idx] = scalar_spectra[idx0];
            idx += nb_scalar_fields;
            all_spectra[idx] = scalar_spectra[idx1];
        }
    }

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
    // call vd2uv to compute u and v in spectral space
    int nb_all_fields = nb_scalar_fields+4*nb_vordiv_fields;
    double all_spectra[2*legendre_size(truncation_)*nb_all_fields];
    for( int jm=0; jm<truncation_; ++jm ) {
        vd2uv(truncation_, jm, nb_vordiv_fields, vorticity_spectra, divergence_spectra, nb_all_fields, all_spectra, nb_scalar_fields, scalar_spectra, config);
    }

    // perform spectral transform to compute all fields in grid point space
    invtrans_uv(nb_all_fields, nb_vordiv_fields, all_spectra, gp_fields, config);

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
