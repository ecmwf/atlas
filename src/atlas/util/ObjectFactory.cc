/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <ostream>

#include "eckit/exception/Exceptions.h"

#include "atlas/runtime/Log.h"
#include "atlas/util/ObjectFactory.h"

using lock_guard = std::lock_guard<std::mutex>;

namespace atlas {
namespace util {

bool ObjectFactoryRegistry::has( const std::string& builder ) const {
    lock_guard lock( mutex_ );
    return ( factories_.find( builder ) != factories_.end() );
}

ObjectFactory* ObjectFactoryRegistry::get( const std::string& builder ) const {
    lock_guard lock( mutex_ );
    auto iterator = factories_.find( builder );

    if ( iterator == factories_.end() ) {
        Log::error() << "No " << factory_ << " for [" << builder << "]" << std::endl;
        Log::error() << "Factories are:" << std::endl;
        for ( const auto& map_pair : factories_ ) {
            Log::error() << "   " << map_pair.first << std::endl;
        }
        throw eckit::SeriousBug( std::string( "No " ) + factory_ + std::string( " called " ) + builder );
    }
    else {
        return iterator->second;
    }
}

void ObjectFactoryRegistry::add( const std::string& builder, ObjectFactory* factory ) {
    lock_guard lock( mutex_ );
    ASSERT( factories_.find( builder ) == factories_.end() );
    factories_[builder] = factory;
    Log::info() << "Registered " << builder << " in " << factory_ << std::endl;
}

void ObjectFactoryRegistry::remove( const std::string& builder ) {
    lock_guard lock( mutex_ );
    factories_.erase( builder );
}

ObjectFactoryRegistry::ObjectFactoryRegistry( const std::string& factory ) : factory_( factory ) {
    Log::info() << "Created Registry" << factory << std::endl;
}

ObjectFactoryRegistry::~ObjectFactoryRegistry() {
    while ( not factories_.empty() ) {
        delete factories_.begin()->second;  // will remove itself from registry in its destructor
    }
}

void ObjectFactoryRegistry::list( std::ostream& out ) const {
    lock_guard lock( mutex_ );
    const char* sep = "";
    for ( const auto& map_pair : factories_ ) {
        out << sep << map_pair.first;
        sep = ", ";
    }
}


//----------------------------------------------------------------------------------------------------------------------

ObjectFactory::ObjectFactory( ObjectFactoryRegistry& registry, const std::string& builder ) :
    registry_( registry ),
    builder_( builder ) {
    if ( not builder_.empty() ) { registry_.add( builder, this ); }
}

ObjectFactory::~ObjectFactory() {
    if ( not builder_.empty() ) { registry_.remove( builder_ ); }
}

//----------------------------------------------------------------------------------------------------------------------

}  // namespace util
}  // namespace atlas
