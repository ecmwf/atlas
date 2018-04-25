/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "eckit/exception/Exceptions.h"
#include "eckit/thread/AutoLock.h"
#include "eckit/thread/Mutex.h"
#include "eckit/utils/Hash.h"

#include "atlas/functionspace.h"
#include "atlas/grid/Grid.h"
#include "atlas/library/defines.h"
#include "atlas/runtime/Log.h"
#include "atlas/trans/Trans.h"

// For factory registration only:
#if ATLAS_HAVE_TRANS
#include "atlas/trans/ifs/TransIFSNodeColumns.h"
#include "atlas/trans/ifs/TransIFSStructuredColumns.h"
#define TRANS_DEFAULT "ifs"
#else
#define TRANS_DEFAULT "local"
#endif
#include "atlas/trans/local_noopt/TransLocal.h"
#include "atlas/trans/localopt/TransLocalopt.h"
#include "atlas/trans/localopt2/TransLocalopt2.h"
#include "atlas/trans/localopt3/TransLocalopt3.h" // --> recommended "local"

namespace atlas {
namespace trans {

TransImpl::~TransImpl() {}

namespace {

static eckit::Mutex* local_mutex               = 0;
static std::map<std::string, TransFactory*>* m = 0;
static pthread_once_t once                     = PTHREAD_ONCE_INIT;

static void init() {
    local_mutex = new eckit::Mutex();
    m           = new std::map<std::string, TransFactory*>();
}

template <typename T>
void load_builder_functionspace() {
    TransBuilderFunctionSpace<T>( "tmp" );
}
template <typename T>
void load_builder_grid() {
    TransBuilderGrid<T>( "tmp" );
}

struct force_link {
    force_link() {
#if ATLAS_HAVE_TRANS
        load_builder_functionspace<TransIFSNodeColumns>();
        load_builder_functionspace<TransIFSStructuredColumns>();
        load_builder_grid<TransIFS>();
#endif
        load_builder_grid<TransLocal>();
        load_builder_grid<TransLocalopt>();
        load_builder_grid<TransLocalopt2>();
        load_builder_grid<TransLocalopt3>();
    }
};

TransFactory& factory( const std::string& name ) {
    std::map<std::string, TransFactory*>::const_iterator j = m->find( name );
    if ( j == m->end() ) {
        Log::error() << "No TransFactory for [" << name << "]" << std::endl;
        Log::error() << "TransFactories are:" << std::endl;
        for ( j = m->begin(); j != m->end(); ++j )
            Log::error() << "   " << ( *j ).first << std::endl;
        throw eckit::SeriousBug( std::string( "No TransFactory called " ) + name );
    }
    return *j->second;
}

}  // namespace

TransFactory::TransFactory( const std::string& name ) : name_( name ) {
    pthread_once( &once, init );

    eckit::AutoLock<eckit::Mutex> lock( local_mutex );

    ASSERT( m->find( name ) == m->end() );
    ( *m )[name] = this;
}

TransFactory::~TransFactory() {
    eckit::AutoLock<eckit::Mutex> lock( local_mutex );
    m->erase( name_ );
}

bool TransFactory::has( const std::string& name ) {
    pthread_once( &once, init );

    eckit::AutoLock<eckit::Mutex> lock( local_mutex );

    static force_link static_linking;

    return ( m->find( name ) != m->end() );
}

void TransFactory::list( std::ostream& out ) {
    pthread_once( &once, init );

    eckit::AutoLock<eckit::Mutex> lock( local_mutex );

    static force_link static_linking;

    const char* sep = "";
    for ( std::map<std::string, TransFactory*>::const_iterator j = m->begin(); j != m->end(); ++j ) {
        out << sep << ( *j ).first;
        sep = ", ";
    }
}

Trans::Implementation* TransFactory::build( const FunctionSpace& gp, const FunctionSpace& sp,
                                            const eckit::Configuration& config ) {
    return build( Cache(), gp, sp, config );
}

Trans::Implementation* TransFactory::build( const Cache& cache, const FunctionSpace& gp, const FunctionSpace& sp,
                                            const eckit::Configuration& config ) {
    pthread_once( &once, init );

    eckit::AutoLock<eckit::Mutex> lock( local_mutex );

    static force_link static_linking;

    std::string suffix( "(" + gp.type() + "," + sp.type() + ")" );
    std::string name = config.getString( "type", TRANS_DEFAULT ) + suffix;

    Log::debug() << "Looking for TransFactory [" << name << "]" << std::endl;

    if ( not config.has( "type" ) and not has( name ) ) {
        name = std::string( "local" ) + suffix;
        Log::debug() << "Looking for TransFactory [" << name << "]" << std::endl;
    }

    return factory( name ).make( cache, gp, sp, config );
}

Trans::Implementation* TransFactory::build( const Grid& grid, int truncation, const eckit::Configuration& config ) {
    return build( Cache(), grid, truncation, config );
}

Trans::Implementation* TransFactory::build( const Cache& cache, const Grid& grid, int truncation,
                                            const eckit::Configuration& config ) {
    pthread_once( &once, init );

    eckit::AutoLock<eckit::Mutex> lock( local_mutex );

    static force_link static_linking;

    std::string name = config.getString( "type", TRANS_DEFAULT );

    Log::debug() << "Looking for TransFactory [" << name << "]" << std::endl;

    if ( not config.has( "type" ) and not has( name ) ) {
        name = std::string( "local" );
        Log::debug() << "Looking for TransFactory [" << name << "]" << std::endl;
    }

    return factory( name ).make( cache, grid, truncation, config );
}

Trans::Trans() {}

Trans::Trans( Implementation* impl ) : impl_( impl ) {}

Trans::Trans( const FunctionSpace& gp, const FunctionSpace& sp, const eckit::Configuration& config ) :
    impl_( TransFactory::build( gp, sp, config ) ) {}

Trans::Trans( const Grid& grid, int truncation, const eckit::Configuration& config ) :
    impl_( TransFactory::build( grid, truncation, config ) ) {}

Trans::Trans( const Cache& cache, const FunctionSpace& gp, const FunctionSpace& sp,
              const eckit::Configuration& config ) :
    impl_( TransFactory::build( cache, gp, sp, config ) ) {}

Trans::Trans( const Cache& cache, const Grid& grid, int truncation, const eckit::Configuration& config ) :
    impl_( TransFactory::build( cache, grid, truncation, config ) ) {}

Trans::Trans( const Trans& trans ) : impl_( trans.impl_ ) {}

int Trans::truncation() const {
    return impl_->truncation();
}

const Grid& Trans::grid() const {
    return impl_->grid();
}

size_t Trans::spectralCoefficients() const {
    return impl_->spectralCoefficients();
}

void Trans::dirtrans( const Field& gpfield, Field& spfield, const eckit::Configuration& config ) const {
    impl_->dirtrans( gpfield, spfield, config );
}

void Trans::dirtrans( const FieldSet& gpfields, FieldSet& spfields, const eckit::Configuration& config ) const {
    impl_->dirtrans( gpfields, spfields, config );
}

void Trans::dirtrans_wind2vordiv( const Field& gpwind, Field& spvor, Field& spdiv,
                                  const eckit::Configuration& config ) const {
    impl_->dirtrans_wind2vordiv( gpwind, spvor, spdiv, config );
}

void Trans::invtrans( const Field& spfield, Field& gpfield, const eckit::Configuration& config ) const {
    impl_->invtrans( spfield, gpfield, config );
}

void Trans::invtrans( const FieldSet& spfields, FieldSet& gpfields, const eckit::Configuration& config ) const {
    impl_->invtrans( spfields, gpfields, config );
}

void Trans::invtrans_grad( const Field& spfield, Field& gradfield, const eckit::Configuration& config ) const {
    impl_->invtrans_grad( spfield, gradfield, config );
}

void Trans::invtrans_grad( const FieldSet& spfields, FieldSet& gradfields, const eckit::Configuration& config ) const {
    impl_->invtrans_grad( spfields, gradfields, config );
}

void Trans::invtrans_vordiv2wind( const Field& spvor, const Field& spdiv, Field& gpwind,
                                  const eckit::Configuration& config ) const {
    impl_->invtrans_vordiv2wind( spvor, spdiv, gpwind, config );
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
    impl_->invtrans( nb_scalar_fields, scalar_spectra, nb_vordiv_fields, vorticity_spectra, divergence_spectra,
                     gp_fields, config );
}

/*!
 * @brief invtrans
 * @param nb_fields
 * @param scalar_spectra
 * @param scalar_fields
 */
void Trans::invtrans( const int nb_scalar_fields, const double scalar_spectra[], double gp_fields[],
                      const eckit::Configuration& config ) const {
    impl_->invtrans( nb_scalar_fields, scalar_spectra, gp_fields, config );
}

/*!
 * @brief Inverse transform of vorticity/divergence to wind(U/V)
 * @param nb_fields [in] Number of fields ( both components of wind count as 1 )
 */
void Trans::invtrans( const int nb_vordiv_fields, const double vorticity_spectra[], const double divergence_spectra[],
                      double gp_fields[], const eckit::Configuration& config ) const {
    impl_->invtrans( nb_vordiv_fields, vorticity_spectra, divergence_spectra, gp_fields, config );
}

/*!
 * @brief Direct transform of scalar fields
 */
void Trans::dirtrans( const int nb_fields, const double scalar_fields[], double scalar_spectra[],
                      const eckit::Configuration& config ) const {
    impl_->dirtrans( nb_fields, scalar_fields, scalar_spectra, config );
}

/*!
 * @brief Direct transform of wind(U/V) to vorticity/divergence
 * @param nb_fields [in] Number of fields ( both components of wind count as 1 )
 */
void Trans::dirtrans( const int nb_fields, const double wind_fields[], double vorticity_spectra[],
                      double divergence_spectra[], const eckit::Configuration& config ) const {
    impl_->dirtrans( nb_fields, wind_fields, vorticity_spectra, divergence_spectra, config );
}

}  // namespace trans
}  // namespace atlas
