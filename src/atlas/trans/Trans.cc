/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "eckit/utils/Hash.h"

#include "atlas/functionspace.h"
#include "atlas/grid/Grid.h"
#include "atlas/library/defines.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"
#include "atlas/trans/Trans.h"
#include "atlas/trans/detail/TransFactory.h"
#include "atlas/trans/detail/TransImpl.h"

namespace atlas {
namespace trans {

namespace {
util::Config options( const eckit::Configuration& config ) {
    util::Config opts = Trans::config();
    opts.set( config );
    return opts;
}
}  // namespace

void Trans::listBackends( std::ostream& out ) {
    TransFactory::list( out );
}

bool Trans::hasBackend( const std::string& backend ) {
    return TransFactory::has( backend );
}

void Trans::backend( const std::string& backend ) {
    ATLAS_ASSERT( hasBackend( backend ) );
    TransFactory::backend( backend );
}

std::string Trans::backend() {
    return TransFactory::backend();
}

const eckit::Configuration& Trans::config() {
    return TransFactory::config();
}

void Trans::config( const eckit::Configuration& _config ) {
    TransFactory::config( _config );
}

Trans::Trans( const FunctionSpace& gp, const FunctionSpace& sp, const eckit::Configuration& config ) :
    Handle( TransFactory::build( gp, sp, config ) ) {}

Trans::Trans( const Grid& grid, int truncation, const eckit::Configuration& config ) :
    Handle( TransFactory::build( grid, truncation, config ) ) {}

Trans::Trans( const Grid& grid, const Domain& domain, int truncation, const eckit::Configuration& config ) :
    Handle( TransFactory::build( grid, domain, truncation, config ) ) {}

Trans::Trans( const Cache& cache, const FunctionSpace& gp, const FunctionSpace& sp,
              const eckit::Configuration& config ) :
    Handle( TransFactory::build( cache, gp, sp, config ) ) {}

Trans::Trans( const Cache& cache, const Grid& grid, int truncation, const eckit::Configuration& config ) :
    Handle( TransFactory::build( cache, grid, truncation, config ) ) {}

Trans::Trans( const Cache& cache, const Grid& grid, const Domain& domain, int truncation,
              const eckit::Configuration& config ) :
    Handle( TransFactory::build( cache, grid, domain, truncation, config ) ) {}

int Trans::truncation() const {
    return get()->truncation();
}

const Grid& Trans::grid() const {
    return get()->grid();
}

const functionspace::Spectral& Trans::spectral() const {
    return get()->spectral();
}

size_t Trans::spectralCoefficients() const {
    return get()->spectralCoefficients();
}

void Trans::dirtrans( const Field& gpfield, Field& spfield, const eckit::Configuration& config ) const {
    get()->dirtrans( gpfield, spfield, options( config ) );
}

void Trans::dirtrans( const FieldSet& gpfields, FieldSet& spfields, const eckit::Configuration& config ) const {
    get()->dirtrans( gpfields, spfields, options( config ) );
}

void Trans::dirtrans_wind2vordiv( const Field& gpwind, Field& spvor, Field& spdiv,
                                  const eckit::Configuration& config ) const {
    get()->dirtrans_wind2vordiv( gpwind, spvor, spdiv, options( config ) );
}

void Trans::invtrans( const Field& spfield, Field& gpfield, const eckit::Configuration& config ) const {
    get()->invtrans( spfield, gpfield, options( config ) );
}

void Trans::invtrans( const FieldSet& spfields, FieldSet& gpfields, const eckit::Configuration& config ) const {
    get()->invtrans( spfields, gpfields, options( config ) );
}

void Trans::invtrans_grad( const Field& spfield, Field& gradfield, const eckit::Configuration& config ) const {
    get()->invtrans_grad( spfield, gradfield, options( config ) );
}

void Trans::invtrans_grad( const FieldSet& spfields, FieldSet& gradfields, const eckit::Configuration& config ) const {
    get()->invtrans_grad( spfields, gradfields, options( config ) );
}

void Trans::invtrans_vordiv2wind( const Field& spvor, const Field& spdiv, Field& gpwind,
                                  const eckit::Configuration& config ) const {
    get()->invtrans_vordiv2wind( spvor, spdiv, gpwind, options( config ) );
}

// -- IFS type fields --
// These fields have special interpretation required. You need to know what
// you're doing.
// See IFS trans library.

/*!
 * @brief invtrans
 * @param nb_scalar_fields
 * @param scalar_spectra
 * @param nb_vordiv_fields
 * @param vorticity_spectra
 * @param divergence_spectra
 * @param gp_fields
 */
void Trans::invtrans( const int nb_scalar_fields, const double scalar_spectra[], const int nb_vordiv_fields,
                      const double vorticity_spectra[], const double divergence_spectra[], double gp_fields[],
                      const eckit::Configuration& config ) const {
    get()->invtrans( nb_scalar_fields, scalar_spectra, nb_vordiv_fields, vorticity_spectra, divergence_spectra,
                     gp_fields, options( config ) );
}

/*!
 * @brief invtrans
 * @param nb_fields
 * @param scalar_spectra
 * @param scalar_fields
 */
void Trans::invtrans( const int nb_scalar_fields, const double scalar_spectra[], double gp_fields[],
                      const eckit::Configuration& config ) const {
    get()->invtrans( nb_scalar_fields, scalar_spectra, gp_fields, options( config ) );
}

/*!
 * @brief Inverse transform of vorticity/divergence to wind(U/V)
 * @param nb_fields [in] Number of fields ( both components of wind count as 1 )
 */
void Trans::invtrans( const int nb_vordiv_fields, const double vorticity_spectra[], const double divergence_spectra[],
                      double gp_fields[], const eckit::Configuration& config ) const {
    get()->invtrans( nb_vordiv_fields, vorticity_spectra, divergence_spectra, gp_fields, options( config ) );
}

/*!
 * @brief Direct transform of scalar fields
 */
void Trans::dirtrans( const int nb_fields, const double scalar_fields[], double scalar_spectra[],
                      const eckit::Configuration& config ) const {
    get()->dirtrans( nb_fields, scalar_fields, scalar_spectra, options( config ) );
}

/*!
 * @brief Direct transform of wind(U/V) to vorticity/divergence
 * @param nb_fields [in] Number of fields ( both components of wind count as 1 )
 */
void Trans::dirtrans( const int nb_fields, const double wind_fields[], double vorticity_spectra[],
                      double divergence_spectra[], const eckit::Configuration& config ) const {
    get()->dirtrans( nb_fields, wind_fields, vorticity_spectra, divergence_spectra, options( config ) );
}

}  // namespace trans
}  // namespace atlas
