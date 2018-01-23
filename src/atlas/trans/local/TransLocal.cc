/*
 * (C) Copyright 2013 ECMWF.
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
