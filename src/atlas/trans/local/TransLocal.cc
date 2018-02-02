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
#include "atlas/trans/VorDivToUV.h"
#include "atlas/option.h"

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

void gp_transpose(
        const int nb_size,
        const int nb_fields,
        const double gp_tmp[],
        double gp_fields[] )
{
    for( int jgp=0; jgp<nb_size; jgp++ ) {
        for( int jfld=0; jfld<nb_fields; jfld++ ) {
            gp_fields[jfld*nb_size+jgp] = gp_tmp[jgp*nb_fields+jfld];
        }
    }
    //Log::info() << "gp_transpose: gp_tmp:" << std::endl;
    //for( int j=0; j<nb_fields*nb_size; j++ ) Log::info() << gp_tmp[j] << " ";
    //Log::info() << std::endl;
    //Log::info() << "gp_transpose: gp_fields:" << std::endl;
    //for( int j=0; j<nb_fields*nb_size; j++ ) Log::info() << gp_fields[j] << " ";
    //Log::info() << std::endl;
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
    std::vector<double> gp_tmp(nb_fields*grid_.size(), 0.);

    // Transform
    if( grid::StructuredGrid g = grid_ ) {
        ATLAS_TRACE( "invtrans_uv structured");
        int idx = 0;
        for( size_t j=0; j<g.ny(); ++j ) {
            double lat = g.y(j) * util::Constants::degreesToRadians();
            double trcFT = fourier_truncation( truncation_, lat );

            // Legendre transform:
            invtrans_legendre( truncation_, trcFT, legPol(lat,j), nb_fields, scalar_spectra, legReal.data(), legImag.data() );

            // Fourier transform:
            for( size_t i=0; i<g.nx(j); ++i ) {
                double lon = g.x(i,j) * util::Constants::degreesToRadians();
                invtrans_fourier( trcFT, lon, nb_fields, legReal.data(), legImag.data(), gp_tmp.data()+(nb_fields*idx));
                for( int jfld=0; jfld<nb_vordiv_fields; ++jfld ) {
                    gp_tmp[nb_fields*idx+jfld] /= std::cos(lat);
                }
                ++idx;
            }
        }
    } else {
        ATLAS_TRACE( "invtrans_uv unstructured");
        int idx = 0;
        for( PointXY p: grid_.xy() ) {
            double lon = p.x() * util::Constants::degreesToRadians();
            double lat = p.y() * util::Constants::degreesToRadians();
            double trcFT = fourier_truncation( truncation_, lat );

            // Legendre transform:
            invtrans_legendre( truncation_, trcFT, legPol(lat,idx), nb_fields, scalar_spectra, legReal.data(), legImag.data() );

            // Fourier transform:
            invtrans_fourier( trcFT, lon, nb_fields, legReal.data(), legImag.data(), gp_tmp.data()+(nb_fields*idx));
            for( int jfld=0; jfld<nb_vordiv_fields; ++jfld ) {
                gp_tmp[nb_fields*idx+jfld] /= std::cos(lat);
            }
            ++idx;
        }
    }

    // transpose result
    gp_transpose( grid_.size(), nb_fields, gp_tmp.data(), gp_fields );

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
    ATLAS_TRACE("TransLocal::invtrans");
    // call vd2uv to compute u and v in spectral space
    int nb_gp = grid_.size();

    int nb_vordiv_spec = 2*legendre_size(truncation_)*nb_vordiv_fields;
    std::vector<double> U(nb_vordiv_spec,0.);
    std::vector<double> V(nb_vordiv_spec,0.);
    trans::VorDivToUV vordiv_to_UV(truncation_, option::type("local"));
    vordiv_to_UV.execute( nb_vordiv_spec, nb_vordiv_fields, vorticity_spectra, divergence_spectra, U.data(), V.data() );

    // perform spectral transform to compute all fields in grid point space
    invtrans_uv(nb_vordiv_fields, nb_vordiv_fields, U.data(), gp_fields, config);
    invtrans_uv(nb_vordiv_fields, nb_vordiv_fields, V.data(), gp_fields+nb_gp*nb_vordiv_fields, config);
    invtrans_uv(nb_scalar_fields, 0, scalar_spectra, gp_fields+2*nb_gp*nb_vordiv_fields, config);

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
