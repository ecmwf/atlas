/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/trans/localopt2/TransLocalopt2.h"
#include "atlas/array.h"
#include "atlas/option.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/runtime/ErrorHandling.h"
#include "atlas/runtime/Log.h"
#include "atlas/trans/VorDivToUV.h"
#include "atlas/trans/local/LegendrePolynomials.h"
#include "atlas/trans/localopt2/FourierTransformsopt2.h"
#include "atlas/trans/localopt2/LegendrePolynomialsopt2.h"
#include "atlas/trans/localopt2/LegendreTransformsopt2.h"
#include "atlas/util/Constants.h"
#include "eckit/linalg/LinearAlgebra.h"
#include "eckit/linalg/Matrix.h"

namespace atlas {
namespace trans {

namespace {
static TransBuilderGrid<TransLocalopt2> builder( "localopt2" );
}

// --------------------------------------------------------------------------------------------------------------------
// Helper functions
// --------------------------------------------------------------------------------------------------------------------
namespace {  // anonymous

size_t legendre_size( const size_t truncation ) {
    return ( truncation + 2 ) * ( truncation + 1 ) / 2;
}

}  // namespace

// --------------------------------------------------------------------------------------------------------------------
// Class TransLocalopt2
// --------------------------------------------------------------------------------------------------------------------

TransLocalopt2::TransLocalopt2( const Cache& cache, const Grid& grid, const long truncation,
                                const eckit::Configuration& config ) :
    grid_( grid ),
    truncation_( truncation ),
    precompute_( config.getBool( "precompute", true ) ) {
    ATLAS_TRACE( "Precompute legendre opt2" );
    int nlats, nlons;
    if ( grid::StructuredGrid( grid_ ) && not grid_.projection() ) {
        grid::StructuredGrid g( grid_ );
        nlats = g.ny();
        nlons = g.nxmax();
    }
    else {
        nlats = grid_.size();
        nlons = grid_.size();
    }
    std::vector<double> lats( nlats );
    std::vector<double> lons( nlons );
    if ( grid::StructuredGrid( grid_ ) && not grid_.projection() ) {
        grid::StructuredGrid g( grid_ );
        // TODO: remove legendre_begin and legendre_data (only legendre_ should be needed)
        for ( size_t j = 0; j < nlats; ++j ) {
            lats[j] = g.y( j ) * util::Constants::degreesToRadians();
        }
        for ( size_t j = 0; j < nlons; ++j ) {
            lons[j] = g.x( j, 0 ) * util::Constants::degreesToRadians();
        }
    }
    else {
        int j( 0 );
        for ( PointXY p : grid_.xy() ) {
            lats[j++] = p.y() * util::Constants::degreesToRadians();
            lons[j++] = p.x() * util::Constants::degreesToRadians();
        }
    }
    // precomputations for Legendre polynomials:
    {
        ATLAS_TRACE( "opt2 precomp Legendre" );
        legendre_.resize( legendre_size( truncation_ + 1 ) * nlats );
        compute_legendre_polynomialsopt2( truncation_ + 1, nlats, lats.data(), legendre_.data() );
    }

    // precomputations for Fourier transformations:
    {
        ATLAS_TRACE( "opt2 precomp Fourier" );
        fourier_.resize( 2 * ( truncation_ + 1 ) * nlons );
        int idx = 0;
        for ( int jlon = 0; jlon < nlons; jlon++ ) {
            for ( int jm = 0; jm < truncation_ + 1; jm++ ) {
                fourier_[idx++] = +std::cos( jm * lons[jlon] );  // real part
                fourier_[idx++] = -std::sin( jm * lons[jlon] );  // imaginary part
            }
        }
    }
    {
        ATLAS_TRACE( "opt2 precomp Fourier tp" );
        fouriertp_.resize( 2 * ( truncation_ + 1 ) * nlons );
        int idx = 0;
        for ( int jm = 0; jm < truncation_ + 1; jm++ ) {
            for ( int jlon = 0; jlon < nlons; jlon++ ) {
                fouriertp_[idx++] = +std::cos( jm * lons[jlon] );  // real part
            }
            for ( int jlon = 0; jlon < nlons; jlon++ ) {
                fouriertp_[idx++] = -std::sin( jm * lons[jlon] );  // imaginary part
            }
        }
    }
}

// --------------------------------------------------------------------------------------------------------------------

TransLocalopt2::TransLocalopt2( const Grid& grid, const long truncation, const eckit::Configuration& config ) :
    TransLocalopt2( Cache(), grid, truncation, config ) {}

// --------------------------------------------------------------------------------------------------------------------

TransLocalopt2::~TransLocalopt2() {}

// --------------------------------------------------------------------------------------------------------------------

void TransLocalopt2::invtrans( const Field& spfield, Field& gpfield, const eckit::Configuration& config ) const {
    NOTIMP;
}

// --------------------------------------------------------------------------------------------------------------------

void TransLocalopt2::invtrans( const FieldSet& spfields, FieldSet& gpfields,
                               const eckit::Configuration& config ) const {
    NOTIMP;
}

// --------------------------------------------------------------------------------------------------------------------

void TransLocalopt2::invtrans_grad( const Field& spfield, Field& gradfield, const eckit::Configuration& config ) const {
    NOTIMP;
}

// --------------------------------------------------------------------------------------------------------------------

void TransLocalopt2::invtrans_grad( const FieldSet& spfields, FieldSet& gradfields,
                                    const eckit::Configuration& config ) const {
    NOTIMP;
}

// --------------------------------------------------------------------------------------------------------------------

void TransLocalopt2::invtrans_vordiv2wind( const Field& spvor, const Field& spdiv, Field& gpwind,
                                           const eckit::Configuration& config ) const {
    NOTIMP;
}

void TransLocalopt2::invtrans( const int nb_scalar_fields, const double scalar_spectra[], double gp_fields[],
                               const eckit::Configuration& config ) const {
    invtrans_uv( truncation_, nb_scalar_fields, 0, scalar_spectra, gp_fields, config );
}

void gp_transposeopt2( const int nb_size, const int nb_fields, const double gp_tmp[], double gp_fields[] ) {
    for ( int jgp = 0; jgp < nb_size; jgp++ ) {
        for ( int jfld = 0; jfld < nb_fields; jfld++ ) {
            gp_fields[jfld * nb_size + jgp] = gp_tmp[jgp * nb_fields + jfld];
        }
    }
}

//-----------------------------------------------------------------------------
// Routine to compute the spectral transform by using a localopt2 Fourier
// transformation
// for a grid (same latitude for all longitudes, allows to compute Legendre
// functions
// once for all longitudes). U and v components are divided by cos(latitude) for
// nb_vordiv_fields > 0.
//
// Author:
// Andreas Mueller *ECMWF*
//
void TransLocalopt2::invtrans_uv( const int truncation, const int nb_scalar_fields, const int nb_vordiv_fields,
                                  const double scalar_spectra[], double gp_fields[],
                                  const eckit::Configuration& config ) const {
    if ( nb_scalar_fields > 0 ) {
        int nb_fields = nb_scalar_fields;

        //eckit::linalg::LinearAlgebra::backend( "string" ) // might want to choose backend with this command

        // Transform
        if ( grid::StructuredGrid g = grid_ ) {
            ATLAS_TRACE( "invtrans_uv structured opt2" );
            int nlats        = g.ny();
            int size_fourier = nb_fields * 2 * g.ny();
            std::vector<double> scl_fourier( size_fourier * ( truncation + 1 ) );

            // Legendre transform:
            {
                ATLAS_TRACE( "opt2 Legendre dgemm" );
                for ( int jm = 0; jm <= truncation; jm++ ) {
                    int noff = ( 2 * truncation + 3 - jm ) * jm / 2, ns = truncation - jm + 1;
                    eckit::linalg::Matrix A( eckit::linalg::Matrix(
                        const_cast<double*>( scalar_spectra ) + nb_fields * 2 * noff, nb_fields * 2, ns ) );
                    eckit::linalg::Matrix B( legendre_.data() + noff * g.ny(), ns, g.ny() );
                    eckit::linalg::Matrix C( scl_fourier.data() + jm * size_fourier, nb_fields * 2, g.ny() );
                    eckit::linalg::LinearAlgebra::backend().gemm( A, B, C );
                }
            }
#if 0  // 1: better for small number of columns, large truncation; 0: better for large number of columns

            // Transposition in Fourier space:
            std::vector<double> scl_fourier_tp( size_fourier * ( truncation + 1 ) );
            {
                ATLAS_TRACE( "opt2 transposition in Fourier" );
                int idx = 0;
                for ( int jm = 0; jm < truncation_ + 1; jm++ ) {
                    for ( int jlat = 0; jlat < g.ny(); jlat++ ) {
                        for ( int imag = 0; imag < 2; imag++ ) {
                            for ( int jfld = 0; jfld < nb_fields; jfld++ ) {
                                int pos_tp = jfld + nb_fields * ( jlat + g.ny() * ( imag + 2 * ( jm ) ) );
                                //int pos  = jfld + nb_fields * ( imag + 2 * ( jlat + g.ny() * ( jm ) ) );
                                scl_fourier_tp[pos_tp] = scl_fourier[idx++];  // = scl_fourier[pos]
                            }
                        }
                    }
                }
            }

            // Fourier transformation:
            std::vector<double> gp_opt2( nb_fields * grid_.size(), 0. );
            {
                ATLAS_TRACE( "opt2 Fourier dgemm" );
                eckit::linalg::Matrix A( scl_fourier_tp.data(), nb_fields * g.ny(), ( truncation_ + 1 ) * 2 );
                eckit::linalg::Matrix B( fourier_.data(), ( truncation_ + 1 ) * 2, g.nxmax() );
                eckit::linalg::Matrix C( gp_opt2.data(), nb_fields * g.ny(), g.nxmax() );
                eckit::linalg::LinearAlgebra::backend().gemm( A, B, C );
            }

            // Transposition in grid point space:
            {
                ATLAS_TRACE( "opt2 transposition in gp-space" );
                int idx = 0;
                for ( int jlon = 0; jlon < g.nxmax(); jlon++ ) {
                    for ( int jlat = 0; jlat < g.ny(); jlat++ ) {
                        for ( int jfld = 0; jfld < nb_fields; jfld++ ) {
                            int pos_tp = jlon + g.nxmax() * ( jlat + g.ny() * ( jfld ) );
                            //int pos  = jfld + nb_fields * ( jlat + g.ny() * ( jlon ) );
                            gp_fields[pos_tp] = gp_opt2[idx++];  // = gp_opt2[pos]
                        }
                    }
                }
            }
#else
            // Transposition in Fourier space:
            std::vector<double> scl_fourier_tp( size_fourier * ( truncation + 1 ) );
            {
                ATLAS_TRACE( "opt2 transposition in Fourier" );
                int idx = 0;
                for ( int jm = 0; jm < truncation_ + 1; jm++ ) {
                    for ( int jlat = 0; jlat < g.ny(); jlat++ ) {
                        for ( int imag = 0; imag < 2; imag++ ) {
                            for ( int jfld = 0; jfld < nb_fields; jfld++ ) {
                                int pos_tp = imag + 2 * ( jm + ( truncation_ + 1 ) * ( jlat + g.ny() * ( jfld ) ) );
                                //int pos  = jfld + nb_fields * ( imag + 2 * ( jlat + g.ny() * ( jm ) ) );
                                scl_fourier_tp[pos_tp] = scl_fourier[idx++];  // = scl_fourier[pos]
                            }
                        }
                    }
                }
            }

            // Fourier transformation:
            std::vector<double> gp_opt2( nb_fields * grid_.size(), 0. );
            {
                ATLAS_TRACE( "opt2 Fourier dgemm" );
                eckit::linalg::Matrix A( fouriertp_.data(), g.nxmax(), ( truncation_ + 1 ) * 2 );
                eckit::linalg::Matrix B( scl_fourier_tp.data(), ( truncation_ + 1 ) * 2, nb_fields * g.ny() );
                eckit::linalg::Matrix C( gp_fields, g.nxmax(), nb_fields * g.ny() );
                eckit::linalg::LinearAlgebra::backend().gemm( A, B, C );
            }

#endif
            // Computing u,v from U,V:
            {
                if ( nb_vordiv_fields > 0 ) {
                    ATLAS_TRACE( "opt2 u,v from U,V" );
                    std::vector<double> coslats( nlats );
                    for ( size_t j = 0; j < nlats; ++j ) {
                        coslats[j] = std::cos( g.y( j ) * util::Constants::degreesToRadians() );
                    }
                    int idx = 0;
                    for ( int jfld = 0; jfld < nb_fields; jfld++ ) {
                        for ( int jlat = 0; jlat < g.ny(); jlat++ ) {
                            for ( int jlon = 0; jlon < g.nxmax(); jlon++ ) {
                                gp_fields[idx] /= coslats[jlat];
                                idx++;
                            }
                        }
                    }
                }
            }
        }
        else {
            ATLAS_TRACE( "invtrans_uv unstructured opt2" );
            int idx = 0;
            for ( PointXY p : grid_.xy() ) {
                double lon   = p.x() * util::Constants::degreesToRadians();
                double lat   = p.y() * util::Constants::degreesToRadians();
                double trcFT = truncation;

                // Legendre transform:
                //invtrans_legendreopt2( truncation, trcFT, truncation_ + 1, legPol( lat, idx ), nb_fields, scalar_spectra,
                //                      legReal.data(), legImag.data() );

                // Fourier transform:
                //invtrans_fourieropt2( trcFT, lon, nb_fields, legReal.data(), legImag.data(),
                //                     gp_tmp.data() + ( nb_fields * idx ) );
                for ( int jfld = 0; jfld < nb_vordiv_fields; ++jfld ) {
                    //gp_tmp[nb_fields * idx + jfld] /= std::cos( lat );
                }
                ++idx;
            }
        }
    }
}

// --------------------------------------------------------------------------------------------------------------------

void TransLocalopt2::invtrans( const int nb_vordiv_fields, const double vorticity_spectra[],
                               const double divergence_spectra[], double gp_fields[],
                               const eckit::Configuration& config ) const {
    invtrans( 0, nullptr, nb_vordiv_fields, vorticity_spectra, divergence_spectra, gp_fields, config );
}

void extend_truncationopt2( const int old_truncation, const int nb_fields, const double old_spectra[],
                            double new_spectra[] ) {
    int k = 0, k_old = 0;
    for ( int m = 0; m <= old_truncation + 1; m++ ) {             // zonal wavenumber
        for ( int n = m; n <= old_truncation + 1; n++ ) {         // total wavenumber
            for ( int imag = 0; imag < 2; imag++ ) {              // imaginary/real part
                for ( int jfld = 0; jfld < nb_fields; jfld++ ) {  // field
                    if ( m == old_truncation + 1 || n == old_truncation + 1 ) { new_spectra[k++] = 0.; }
                    else {
                        new_spectra[k++] = old_spectra[k_old++];
                    }
                }
            }
        }
    }
}

// --------------------------------------------------------------------------------------------------------------------

void TransLocalopt2::invtrans( const int nb_scalar_fields, const double scalar_spectra[], const int nb_vordiv_fields,
                               const double vorticity_spectra[], const double divergence_spectra[], double gp_fields[],
                               const eckit::Configuration& config ) const {
    ATLAS_TRACE( "TransLocalopt2::invtrans" );
    int nb_gp              = grid_.size();
    int nb_vordiv_spec_ext = 2 * legendre_size( truncation_ + 1 ) * nb_vordiv_fields;
    if ( nb_vordiv_fields > 0 ) {
        std::vector<double> vorticity_spectra_extended( nb_vordiv_spec_ext, 0. );
        std::vector<double> divergence_spectra_extended( nb_vordiv_spec_ext, 0. );
        std::vector<double> U_ext( nb_vordiv_spec_ext, 0. );
        std::vector<double> V_ext( nb_vordiv_spec_ext, 0. );

        {
            ATLAS_TRACE( "opt2 extend vordiv" );
            // increase truncation in vorticity_spectra and divergence_spectra:
            extend_truncationopt2( truncation_, nb_vordiv_fields, vorticity_spectra,
                                   vorticity_spectra_extended.data() );
            extend_truncationopt2( truncation_, nb_vordiv_fields, divergence_spectra,
                                   divergence_spectra_extended.data() );
        }

        {
            ATLAS_TRACE( "vordiv to UV opt2" );
            // call vd2uv to compute u and v in spectral space
            trans::VorDivToUV vordiv_to_UV_ext( truncation_ + 1, option::type( "localopt2" ) );
            vordiv_to_UV_ext.execute( nb_vordiv_spec_ext, nb_vordiv_fields, vorticity_spectra_extended.data(),
                                      divergence_spectra_extended.data(), U_ext.data(), V_ext.data() );
        }

        // perform spectral transform to compute all fields in grid point space
        invtrans_uv( truncation_ + 1, nb_vordiv_fields, nb_vordiv_fields, U_ext.data(), gp_fields, config );
        invtrans_uv( truncation_ + 1, nb_vordiv_fields, nb_vordiv_fields, V_ext.data(),
                     gp_fields + nb_gp * nb_vordiv_fields, config );
    }
    if ( nb_scalar_fields > 0 ) {
        int nb_scalar_spec_ext = 2 * legendre_size( truncation_ + 1 ) * nb_scalar_fields;
        std::vector<double> scalar_spectra_extended( nb_scalar_spec_ext, 0. );
        extend_truncationopt2( truncation_, nb_scalar_fields, scalar_spectra, scalar_spectra_extended.data() );
        invtrans_uv( truncation_ + 1, nb_scalar_fields, 0, scalar_spectra_extended.data(),
                     gp_fields + 2 * nb_gp * nb_vordiv_fields, config );
    }
}

// --------------------------------------------------------------------------------------------------------------------

void TransLocalopt2::dirtrans( const Field& gpfield, Field& spfield, const eckit::Configuration& config ) const {
    NOTIMP;
    // Not implemented and not planned.
    // Use the TransIFS implementation instead.
}

// --------------------------------------------------------------------------------------------------------------------

void TransLocalopt2::dirtrans( const FieldSet& gpfields, FieldSet& spfields,
                               const eckit::Configuration& config ) const {
    NOTIMP;
    // Not implemented and not planned.
    // Use the TransIFS implementation instead.
}

// --------------------------------------------------------------------------------------------------------------------

void TransLocalopt2::dirtrans_wind2vordiv( const Field& gpwind, Field& spvor, Field& spdiv,
                                           const eckit::Configuration& config ) const {
    NOTIMP;
    // Not implemented and not planned.
    // Use the TransIFS implementation instead.
}

// --------------------------------------------------------------------------------------------------------------------

void TransLocalopt2::dirtrans( const int nb_fields, const double scalar_fields[], double scalar_spectra[],
                               const eckit::Configuration& ) const {
    NOTIMP;
    // Not implemented and not planned.
    // Use the TransIFS implementation instead.
}

// --------------------------------------------------------------------------------------------------------------------

void TransLocalopt2::dirtrans( const int nb_fields, const double wind_fields[], double vorticity_spectra[],
                               double divergence_spectra[], const eckit::Configuration& ) const {
    NOTIMP;
    // Not implemented and not planned.
    // Use the TransIFS implementation instead.
}

// --------------------------------------------------------------------------------------------------------------------

}  // namespace trans
}  // namespace atlas
