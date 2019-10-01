/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/field/State.h"

#include <iomanip>
#include <map>
#include <memory>
#include <string>

#include "eckit/thread/AutoLock.h"
#include "eckit/thread/Mutex.h"

#include "atlas/field/Field.h"
#include "atlas/grid/Grid.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"

namespace atlas {
namespace field {

namespace {

static eckit::Mutex* local_mutex                        = nullptr;
static std::map<std::string, StateGeneratorFactory*>* m = nullptr;
static pthread_once_t once                              = PTHREAD_ONCE_INIT;

static void init() {
    local_mutex = new eckit::Mutex();
    m           = new std::map<std::string, StateGeneratorFactory*>();
}

template <typename T>
void load_builder() {
    StateGeneratorBuilder<T>( "tmp" );
}

struct force_link {
    force_link() = default;
};
}  // namespace

void State::initialize( const std::string& generator, const eckit::Parametrisation& params ) {
    std::unique_ptr<StateGenerator> state_generator( StateGeneratorFactory::build( generator, params ) );
    state_generator->generate( *this, params );
}

//------------------------------------------------------------------------------------------------------

State::State() = default;

State::State( const std::string& generator, const eckit::Parametrisation& params ) {
    initialize( generator, params );
}

const util::Metadata& State::metadata() const {
    return metadata_;
}

util::Metadata& State::metadata() {
    return metadata_;
}

Field State::add( Field field ) {
    ATLAS_ASSERT( field );

    if ( field.name().empty() ) {
        std::stringstream new_name;
        new_name << "field_" << std::setw( 5 ) << std::setfill( '0' ) << fields_.size();
        ATLAS_ASSERT( !has( new_name.str() ) );
        field.rename( new_name.str() );
    }

    if ( has( field.name() ) ) {
        std::stringstream msg;
        msg << "Trying to add field '" << field.name() << "' to State, but State already has a field with this name.";
        throw_AssertionFailed( msg.str(), Here() );
    }
    fields_[field.name()] = field;
    return field;
}

const Field& State::field( const std::string& name ) const {
    if ( !has( name ) ) {
        std::ostringstream msg;
        msg << "Trying to access field `" << name << "' in State, but no field with this name is present in State.";
        throw_AssertionFailed( msg.str(), Here() );
    }
    return fields_.find( name )->second;
}

Field& State::field( const std::string& name ) {
    return const_cast<Field&>( static_cast<const State*>( this )->field( name ) );
}

const Field& State::field( const idx_t idx ) const {
    if ( idx >= size() ) {
        std::ostringstream msg;
        msg << "Trying to access field in State with index " << idx << ", but there exist only " << size()
            << " fields in State.";
        throw_AssertionFailed( msg.str(), Here() );
    }
    FieldMap::const_iterator it = fields_.begin();
    for ( idx_t i = 0; i < idx; ++i ) {
        ++it;
    }
    return it->second;
}

Field& State::field( const idx_t idx ) {
    return const_cast<Field&>( static_cast<const State*>( this )->field( idx ) );
}

std::vector<std::string> State::field_names() const {
    std::vector<std::string> ret;
    if ( fields_.size() ) {
        ret.reserve( fields_.size() );
    }

    for ( FieldMap::const_iterator it = fields_.begin(); it != fields_.end(); ++it ) {
        ret.push_back( it->first );
    }
    return ret;
}

void State::remove( const std::string& name ) {
    if ( fields_.find( name ) == fields_.end() ) {
        std::stringstream msg;
        msg << "Trying to remove field '" << name << "' from State, but it is not present in State.";
        throw_AssertionFailed( msg.str(), Here() );
    }
    fields_.erase( name );
}

//-----------------------------------------------------------------------------

StateGenerator::StateGenerator( const eckit::Parametrisation& ) {}

StateGenerator::~StateGenerator() = default;

StateGenerator* StateGeneratorFactory::build( const std::string& name, const eckit::Parametrisation& param ) {
    pthread_once( &once, init );

    eckit::AutoLock<eckit::Mutex> lock( local_mutex );

    static force_link static_linking;

    std::map<std::string, StateGeneratorFactory*>::const_iterator j = m->find( name );

    Log::debug() << "Looking for StateGeneratorFactory [" << name << "]" << std::endl;

    if ( j == m->end() ) {
        Log::error() << "No StateGeneratorFactory for [" << name << "]" << std::endl;
        Log::error() << "StateFactories are:" << std::endl;
        for ( j = m->begin(); j != m->end(); ++j ) {
            Log::error() << "   " << ( *j ).first << std::endl;
        }
        throw_Exception( std::string( "No StateGeneratorFactory called " ) + name, Here() );
    }

    return ( *j ).second->make( param );
}

void StateGeneratorFactory::list( std::ostream& out ) {
    pthread_once( &once, init );

    eckit::AutoLock<eckit::Mutex> lock( local_mutex );

    static force_link static_linking;

    const char* sep = "";
    for ( std::map<std::string, StateGeneratorFactory*>::const_iterator j = m->begin(); j != m->end(); ++j ) {
        out << sep << ( *j ).first;
        sep = ", ";
    }
}

bool StateGeneratorFactory::has( const std::string& name ) {
    pthread_once( &once, init );

    eckit::AutoLock<eckit::Mutex> lock( local_mutex );

    static force_link static_linking;

    return ( m->find( name ) != m->end() );
}

StateGeneratorFactory::StateGeneratorFactory( const std::string& name ) : name_( name ) {
    pthread_once( &once, init );

    eckit::AutoLock<eckit::Mutex> lock( local_mutex );

    ATLAS_ASSERT( m->find( name ) == m->end() );
    ( *m )[name] = this;
}

StateGeneratorFactory::~StateGeneratorFactory() {
    eckit::AutoLock<eckit::Mutex> lock( local_mutex );
    m->erase( name_ );
}

//-----------------------------------------------------------------------------
// C wrapper interfaces to C++ routines

extern "C" {

State* atlas__State__new() {
    return new State;
}

void atlas__State__initialize( State* This, const char* generator, const eckit::Parametrisation* params ) {
    ATLAS_ASSERT( This != nullptr, "Reason: Use of uninitialised atlas_State" );
    This->initialize( std::string( generator ), *params );
}

void atlas__State__delete( State* This ) {
    ATLAS_ASSERT( This != nullptr, "Reason: Use of uninitialised atlas_State" );
    delete This;
}

void atlas__State__add( State* This, FieldImpl* field ) {
    ATLAS_ASSERT( This != nullptr, "Reason: Use of uninitialised atlas_State" );
    This->add( field );
}

void atlas__State__remove( State* This, const char* name ) {
    ATLAS_ASSERT( This != nullptr, "Reason: Use of uninitialised atlas_State" );
    This->remove( name );
}

int atlas__State__has( State* This, const char* name ) {
    ATLAS_ASSERT( This != nullptr, "Reason: Use of uninitialised atlas_State" );
    return This->has( name );
}

FieldImpl* atlas__State__field_by_name( State* This, const char* name ) {
    ATLAS_ASSERT( This != nullptr, "Reason: Use of uninitialised atlas_State" );
    return This->field( std::string( name ) ).get();
}

FieldImpl* atlas__State__field_by_index( State* This, idx_t index ) {
    ATLAS_ASSERT( This != nullptr, "Reason: Use of uninitialised atlas_State" );
    return This->field( index ).get();
}

idx_t atlas__State__size( const State* This ) {
    ATLAS_ASSERT( This != nullptr, "Reason: Use of uninitialised atlas_State" );
    return This->size();
}

util::Metadata* atlas__State__metadata( State* This ) {
    ATLAS_ASSERT( This != nullptr, "Reason: Use of uninitialised atlas_State" );
    return &This->metadata();
}
}
//-----------------------------------------------------------------------------

}  // namespace field
}  // namespace atlas
