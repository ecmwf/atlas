/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/grid/Spec.h"

#include <iostream>
#include <mutex>

#include "atlas/grid/Grid.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"
#include "atlas/util/Config.h"

namespace atlas {

Spec::Spec( const std::string& id ) : Handle( nullptr /*SpecFactory::create( id )*/ ) {}

std::string Spec::name() const {
    ATLAS_ASSERT_MSG( operator bool(), "Spec::name(): not assigned" );
    return get()->getString( "name" );
}

std::string Spec::type() const {
    ATLAS_ASSERT_MSG( operator bool(), "Spec::type(): not assigned" );
    return get()->getString( "type" );
}

std::string Spec::uid() const {
    ATLAS_ASSERT_MSG( operator bool(), "Spec::uid(): not assigned" );
    return get()->has( "uid" ) ? get()->getString( "uid" ) : Grid( *get() ).uid();
}

void Spec::hash( eckit::Hash& h ) const {
    ATLAS_ASSERT_MSG( operator bool(), "Spec::hash(): not assigned" );
    Grid( *get() ).hash( h );
}

//---------------------------------------------------------------------------------------------------------------------

static std::once_flag once;
static std::recursive_mutex* mtx                            = nullptr;
static std::map<std::string, const Spec::Implementation>* m = nullptr;

static void init() {
    mtx = new std::recursive_mutex();
    m   = new std::map<std::string, const Spec::Implementation>();
}

SpecFactory::SpecFactory( const std::string& id, const Spec::Implementation& spec ) {
    std::call_once( once, init );
    std::lock_guard<std::recursive_mutex> lock( *mtx );

    if ( m->find( id ) != m->end() ) {
        throw_Exception( "SpecFactory: duplicate '" + id + "'" );
    }

    m->emplace( id, spec );
}

SpecFactory::SpecFactory( const std::string& id, const eckit::PathName& path ) :
    SpecFactory( id, Spec::Implementation( path ) ) {}

Spec::Implementation* SpecFactory::create( const std::string& id ) {
    std::call_once( once, init );
    std::lock_guard<std::recursive_mutex> lock( *mtx );

    auto j = m->find( id );
    if ( j != m->end() ) {
        return new Spec::Implementation( j->second );
    }

    list( Log::error() << "SpecFactory: unknown '" << id << "', choices are: " );
    throw_Exception( "SpecFactory: unknown '" + id + "'" );
}

void SpecFactory::list( std::ostream& out ) {
    std::call_once( once, init );
    std::lock_guard<std::recursive_mutex> lock( *mtx );

    const char* sep = "";
    for ( const auto& j : *m ) {
        out << sep << j.first;
        sep = ", ";
    }
}

}  // namespace atlas
