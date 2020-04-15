/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <map>
#include <string>

#include "eckit/thread/AutoLock.h"
#include "eckit/thread/Mutex.h"

#include "atlas/library/config.h"
#include "atlas/numerics/Method.h"
#include "atlas/numerics/Nabla.h"
#include "atlas/numerics/fvm/Method.h"
#include "atlas/numerics/fvm/Nabla.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"
#include "atlas/util/Config.h"

namespace {

static eckit::Mutex* local_mutex                                = nullptr;
static std::map<std::string, atlas::numerics::NablaFactory*>* m = nullptr;
static pthread_once_t once                                      = PTHREAD_ONCE_INIT;

static void init() {
    local_mutex = new eckit::Mutex();
    m           = new std::map<std::string, atlas::numerics::NablaFactory*>();
}
}  // namespace

namespace atlas {
namespace numerics {

NablaImpl::NablaImpl( const Method& method, const eckit::Parametrisation& ) : method_( &method ) {}

NablaImpl::~NablaImpl() = default;

Nabla::Nabla( const Method& method, const eckit::Parametrisation& p ) : Handle( NablaFactory::build( method, p ) ) {}

Nabla::Nabla( const Method& method ) : Nabla( method, util::NoConfig() ) {}

void Nabla::gradient( const Field& scalar, Field& grad ) const {
    get()->gradient( scalar, grad );
}

void Nabla::divergence( const Field& vector, Field& div ) const {
    get()->divergence( vector, div );
}

void Nabla::curl( const Field& vector, Field& curl ) const {
    get()->curl( vector, curl );
}

void Nabla::laplacian( const Field& scalar, Field& laplacian ) const {
    get()->laplacian( scalar, laplacian );
}

namespace {

template <typename T>
void load_builder() {
    NablaBuilder<T>( "tmp" );
}

struct force_link {
    force_link() { load_builder<fvm::Nabla>(); }
};

}  // namespace

NablaFactory::NablaFactory( const std::string& name ) : name_( name ) {
    pthread_once( &once, init );

    eckit::AutoLock<eckit::Mutex> lock( local_mutex );

    ATLAS_ASSERT( m->find( name ) == m->end() );
    ( *m )[name] = this;
}

NablaFactory::~NablaFactory() {
    eckit::AutoLock<eckit::Mutex> lock( local_mutex );
    m->erase( name_ );
}

void NablaFactory::list( std::ostream& out ) {
    pthread_once( &once, init );

    eckit::AutoLock<eckit::Mutex> lock( local_mutex );

    static force_link static_linking;

    const char* sep = "";
    for ( std::map<std::string, NablaFactory*>::const_iterator j = m->begin(); j != m->end(); ++j ) {
        out << sep << ( *j ).first;
        sep = ", ";
    }
}

bool NablaFactory::has( const std::string& name ) {
    pthread_once( &once, init );

    eckit::AutoLock<eckit::Mutex> lock( local_mutex );

    static force_link static_linking;

    return ( m->find( name ) != m->end() );
}

const NablaImpl* NablaFactory::build( const Method& method, const eckit::Parametrisation& p ) {
    pthread_once( &once, init );

    eckit::AutoLock<eckit::Mutex> lock( local_mutex );

    static force_link static_linking;

    std::map<std::string, NablaFactory*>::const_iterator j = m->find( method.name() );

    Log::debug() << "Looking for NablaFactory [" << method.name() << "]" << '\n';

    if ( j == m->end() ) {
        Log::error() << "No NablaFactory for [" << method.name() << "]" << '\n';
        Log::error() << "NablaFactories are:" << '\n';
        for ( j = m->begin(); j != m->end(); ++j ) {
            Log::error() << "   " << ( *j ).first << '\n';
        }
        throw_Exception( std::string( "No NablaFactory called " ) + method.name() );
    }

    return ( *j ).second->make( method, p );
}

extern "C" {

void atlas__Nabla__delete( Nabla::Implementation* This ) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialisd atlas_numerics_Nabla" );
    delete This;
}

const Nabla::Implementation* atlas__Nabla__create( const Method* method, const eckit::Parametrisation* config ) {
    ATLAS_ASSERT( method != nullptr, "Cannot access uninitialisd atlas_numerics_Method" );
    ATLAS_ASSERT( config != nullptr, "Cannot access uninitialisd atlas_Config" );
    const Nabla::Implementation* nabla( nullptr );
    {
        Nabla n( *method, *config );
        nabla = n.get();
        nabla->attach();
    }
    nabla->detach();
    return nabla;
}

void atlas__Nabla__gradient( const Nabla::Implementation* This, const field::FieldImpl* scalar,
                             field::FieldImpl* grad ) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialisd atlas_numerics_Nabla" );
    ATLAS_ASSERT( scalar != nullptr, "Cannot access uninitialisd atlas_Field" );
    ATLAS_ASSERT( grad != nullptr, "Cannot access uninitialisd atlas_Field" );
    Field fgrad( grad );
    This->gradient( scalar, fgrad );
}

void atlas__Nabla__divergence( const Nabla::Implementation* This, const field::FieldImpl* vector,
                               field::FieldImpl* div ) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialisd atlas_numerics_Nabla" );
    ATLAS_ASSERT( vector != nullptr, "Cannot access uninitialisd atlas_Field" );
    ATLAS_ASSERT( div != nullptr, "Cannot access uninitialisd atlas_Field" );
    Field fdiv( div );
    This->divergence( vector, fdiv );
}

void atlas__Nabla__curl( const Nabla::Implementation* This, const field::FieldImpl* vector, field::FieldImpl* curl ) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialisd atlas_numerics_Nabla" );
    ATLAS_ASSERT( vector != nullptr, "Cannot access uninitialisd atlas_Field" );
    ATLAS_ASSERT( curl != nullptr, "Cannot access uninitialisd atlas_Field" );
    Field fcurl( curl );
    This->curl( vector, fcurl );
}

void atlas__Nabla__laplacian( const Nabla::Implementation* This, const field::FieldImpl* scalar,
                              field::FieldImpl* laplacian ) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialisd atlas_numerics_Nabla" );
    ATLAS_ASSERT( scalar != nullptr, "Cannot access uninitialisd atlas_Field" );
    ATLAS_ASSERT( laplacian != nullptr, "Cannot access uninitialisd atlas_Field" );
    Field flaplacian( laplacian );
    This->laplacian( scalar, flaplacian );
}

const functionspace::FunctionSpaceImpl* atlas__Nabla__functionspace( const Nabla::Implementation* This ) {
    return This->functionspace().get();
}

}  // extern "C"

}  // namespace numerics
}  // namespace atlas
