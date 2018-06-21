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

namespace {
struct default_backend {
#if ATLAS_HAVE_TRANS
    std::string value = "ifs";
#else
    std::string value = "local";
#endif
    static default_backend instance() {
        static default_backend x;
        return x;
    }

private:
    default_backend() = default;
};
}  // namespace

namespace atlas {
namespace trans {


class TransBackend {
public:
    static void list( std::ostream& );
    static bool has( const std::string& name );
    static void backend( const std::string& );
    static std::string backend();
    static void config( const eckit::Configuration& );
    static const eckit::Configuration& config();
private:
    static util::Config default_options_;
};

util::Config TransBackend::default_options_ = util::Config( "type", default_backend::instance().value );

TransImpl::~TransImpl() {}

namespace {

static eckit::Mutex* local_mutex               = 0;
static std::map<std::string, int> *b = 0;
static std::map<std::string, TransFactory*>* m = 0;
static pthread_once_t once                     = PTHREAD_ONCE_INIT;

static void init() {
    local_mutex = new eckit::Mutex();
    b           = new std::map<std::string, int>();
    m           = new std::map<std::string, TransFactory*>();
}

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

TransFactory::TransFactory( const std::string& name, const std::string& backend ) : name_( name ), backend_( backend ) {
    pthread_once( &once, init );

    eckit::AutoLock<eckit::Mutex> lock( local_mutex );

    ASSERT( m->find( name ) == m->end() );
    ( *m )[name] = this;

    if( b->find( backend ) == b->end() ) ( *b )[backend] = 0;
    ( *b )[backend]++;
}

TransFactory::~TransFactory() {
    eckit::AutoLock<eckit::Mutex> lock( local_mutex );
    m->erase( name_ );
    ( *b )[backend_]--;
    if ( ( *b )[backend_] == 0 ) b->erase( backend_ );
}

bool TransFactory::has( const std::string& name ) {
    pthread_once( &once, init );
    eckit::AutoLock<eckit::Mutex> lock( local_mutex );
    return ( m->find( name ) != m->end() );
}

bool TransBackend::has( const std::string& name ) {
    pthread_once( &once, init );
    eckit::AutoLock<eckit::Mutex> lock( local_mutex );
    return ( b->find( name ) != b->end() );
}

void TransBackend::backend( const std::string& backend ) {
    pthread_once( &once, init );
    eckit::AutoLock<eckit::Mutex> lock( local_mutex );
    default_options_.set( "type", backend );
}

std::string TransBackend::backend() {
    return default_options_.getString( "type" );
}

const eckit::Configuration& TransBackend::config() {
    return default_options_;
}

void TransBackend::config( const eckit::Configuration& config ) {
    std::string type = default_options_.getString( "type" );
    default_options_ = config;
    if ( not config.has( "type" ) ) { default_options_.set( "type", type ); }
}

void TransBackend::list( std::ostream& out ) {
    pthread_once( &once, init );
    eckit::AutoLock<eckit::Mutex> lock( local_mutex );
    const char* sep = "";
    for( std::map<std::string,int>::const_iterator j = b->begin(); j != b->end(); ++j ) {
        out << sep << ( *j ).first;
        sep = ", ";
    }
}

void TransFactory::list( std::ostream& out ) {
    TransBackend::list( out );
}

Trans::Implementation* TransFactory::build( const FunctionSpace& gp, const FunctionSpace& sp,
                                            const eckit::Configuration& config ) {
    return build( Cache(), gp, sp, config );
}

Trans::Implementation* TransFactory::build( const Cache& cache, const FunctionSpace& gp, const FunctionSpace& sp,
                                            const eckit::Configuration& config ) {
    if ( cache.trans() ) {
        Log::debug() << "Creating Trans from cache, ignoring any other arguments" << std::endl;
        return cache.trans();
    }

    pthread_once( &once, init );

    eckit::AutoLock<eckit::Mutex> lock( local_mutex );

    util::Config options = TransBackend::config();
    options.set( config );

    std::string backend = options.getString( "type" );

    Log::debug() << "Looking for TransFactory [" << backend << "]" << std::endl;
    if( ! TransBackend::has( backend ) ) {
      Log::error() << "No TransFactory for [" << backend << "]" << std::endl;
      Log::error() << "TransFactories are :" << std::endl;
      TransBackend::list( Log::error() );
      Log::error() << std::endl;
      throw eckit::SeriousBug( std::string( "No TransFactory called ") + backend );
    }

    std::string suffix( "(" + gp.type() + "," + sp.type() + ")" );
    std::string name = backend + suffix;

    Log::debug() << "Looking for TransFactory [" << name << "]" << std::endl;

    return factory( name ).make( cache, gp, sp, options );
}

Trans::Implementation* TransFactory::build( const Grid& grid, int truncation, const eckit::Configuration& config ) {
    return build( Cache(), grid, truncation, config );
}

Trans::Implementation* TransFactory::build( const Grid& grid, const Domain& domain, int truncation,
                                            const eckit::Configuration& config ) {
    return build( Cache(), grid, domain, truncation, config );
}

Trans::Implementation* TransFactory::build( const Cache& cache, const Grid& grid, int truncation,
                                            const eckit::Configuration& config ) {
    return build( cache, grid, grid.domain(), truncation, config );
}

Trans::Implementation* TransFactory::build( const Cache& cache, const Grid& grid, const Domain& domain, int truncation,
                                            const eckit::Configuration& config ) {
    if ( cache.trans() ) {
        Log::debug() << "Creating Trans from cache, ignoring any other arguments" << std::endl;
        return cache.trans();
    }

    pthread_once( &once, init );

    eckit::AutoLock<eckit::Mutex> lock( local_mutex );

    util::Config options = TransBackend::config();
    options.set( config );

    std::string backend = options.getString( "type" );

    Log::debug() << "Looking for TransFactory [" << backend << "]" << std::endl;
    if( ! TransBackend::has( backend ) ) {
        Log::error() << "No TransFactory for [" << backend << "]" << std::endl;
        Log::error() << "TransFactories are :" << std::endl;
        TransBackend::list( Log::error() );
        Log::error() << std::endl;
        throw eckit::SeriousBug( std::string( "No TransFactory called ") + backend );
    }

    std::string name = backend;

    return factory( name ).make( cache, grid, domain, truncation, options );
}

bool Trans::hasBackend( const std::string& backend ) {
    return TransBackend::has( backend );
}

void Trans::backend( const std::string& backend ) {
    ASSERT( hasBackend( backend ) );
    TransBackend::backend( backend );
}

std::string Trans::backend() {
    return TransBackend::backend();
}

const eckit::Configuration& Trans::config() {
    return TransBackend::config();
}

void Trans::config( const eckit::Configuration& options ) {
    TransBackend::config( options );
}

namespace {
util::Config options( const eckit::Configuration& config ) {
    util::Config opts = Trans::config();
    opts.set( config );
    return opts;
}
}  // namespace

Trans::Trans() {}

Trans::Trans( Implementation* impl ) : impl_( impl ) {}

Trans::Trans( const FunctionSpace& gp, const FunctionSpace& sp, const eckit::Configuration& config ) :
    impl_( TransFactory::build( gp, sp, config ) ) {}

Trans::Trans( const Grid& grid, int truncation, const eckit::Configuration& config ) :
    impl_( TransFactory::build( grid, truncation, config ) ) {}

Trans::Trans( const Grid& grid, const Domain& domain, int truncation, const eckit::Configuration& config ) :
    impl_( TransFactory::build( grid, domain, truncation, config ) ) {}

Trans::Trans( const Cache& cache, const FunctionSpace& gp, const FunctionSpace& sp,
              const eckit::Configuration& config ) :
    impl_( TransFactory::build( cache, gp, sp, config ) ) {}

Trans::Trans( const Cache& cache, const Grid& grid, int truncation, const eckit::Configuration& config ) :
    impl_( TransFactory::build( cache, grid, truncation, config ) ) {}

Trans::Trans( const Cache& cache, const Grid& grid, const Domain& domain, int truncation,
              const eckit::Configuration& config ) :
    impl_( TransFactory::build( cache, grid, domain, truncation, config ) ) {}

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
    impl_->dirtrans( gpfield, spfield, options( config ) );
}

void Trans::dirtrans( const FieldSet& gpfields, FieldSet& spfields, const eckit::Configuration& config ) const {
    impl_->dirtrans( gpfields, spfields, options( config ) );
}

void Trans::dirtrans_wind2vordiv( const Field& gpwind, Field& spvor, Field& spdiv,
                                  const eckit::Configuration& config ) const {
    impl_->dirtrans_wind2vordiv( gpwind, spvor, spdiv, options( config ) );
}

void Trans::invtrans( const Field& spfield, Field& gpfield, const eckit::Configuration& config ) const {
    impl_->invtrans( spfield, gpfield, options( config ) );
}

void Trans::invtrans( const FieldSet& spfields, FieldSet& gpfields, const eckit::Configuration& config ) const {
    impl_->invtrans( spfields, gpfields, options( config ) );
}

void Trans::invtrans_grad( const Field& spfield, Field& gradfield, const eckit::Configuration& config ) const {
    impl_->invtrans_grad( spfield, gradfield, options( config ) );
}

void Trans::invtrans_grad( const FieldSet& spfields, FieldSet& gradfields, const eckit::Configuration& config ) const {
    impl_->invtrans_grad( spfields, gradfields, options( config ) );
}

void Trans::invtrans_vordiv2wind( const Field& spvor, const Field& spdiv, Field& gpwind,
                                  const eckit::Configuration& config ) const {
    impl_->invtrans_vordiv2wind( spvor, spdiv, gpwind, options( config ) );
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
    impl_->invtrans( nb_scalar_fields, scalar_spectra, gp_fields, options( config ) );
}

/*!
 * @brief Inverse transform of vorticity/divergence to wind(U/V)
 * @param nb_fields [in] Number of fields ( both components of wind count as 1 )
 */
void Trans::invtrans( const int nb_vordiv_fields, const double vorticity_spectra[], const double divergence_spectra[],
                      double gp_fields[], const eckit::Configuration& config ) const {
    impl_->invtrans( nb_vordiv_fields, vorticity_spectra, divergence_spectra, gp_fields, options( config ) );
}

/*!
 * @brief Direct transform of scalar fields
 */
void Trans::dirtrans( const int nb_fields, const double scalar_fields[], double scalar_spectra[],
                      const eckit::Configuration& config ) const {
    impl_->dirtrans( nb_fields, scalar_fields, scalar_spectra, options( config ) );
}

/*!
 * @brief Direct transform of wind(U/V) to vorticity/divergence
 * @param nb_fields [in] Number of fields ( both components of wind count as 1 )
 */
void Trans::dirtrans( const int nb_fields, const double wind_fields[], double vorticity_spectra[],
                      double divergence_spectra[], const eckit::Configuration& config ) const {
    impl_->dirtrans( nb_fields, wind_fields, vorticity_spectra, divergence_spectra, options( config ) );
}

}  // namespace trans
}  // namespace atlas
