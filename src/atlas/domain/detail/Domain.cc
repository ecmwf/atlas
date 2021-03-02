/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <string>

#include "eckit/utils/MD5.h"

#include "atlas/domain/detail/Domain.h"
#include "atlas/domain/detail/DomainFactory.h"
#include "atlas/projection/Projection.h"
#include "atlas/runtime/Exception.h"

namespace atlas {
namespace domain {

const Domain* Domain::create() {
    // default: global domain
    util::Config projParams;
    projParams.set( "type", "global" );
    return Domain::create( projParams );
}

const Domain* Domain::create( const eckit::Parametrisation& p ) {
    std::string domain_type;
    if ( p.get( "type", domain_type ) ) {
        return DomainFactory::build( domain_type, p );
    }

    // should return error here
    throw_Exception( "type missing in Params", Here() );
}

//---------------------------------------------------------------------------------------------------------------------

extern "C" {
const Domain* atlas__Domain__ctor_config( const eckit::Parametrisation* config ) {
    return Domain::create( *config );
}
void atlas__Domain__type( const Domain* This, char*& type, int& size ) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialised atlas_Domain" );
    std::string s = This->type();
    size          = static_cast<int>( s.size() + 1 );
    type          = new char[size];
    strcpy( type, s.c_str() );
}
void atlas__Domain__hash( const Domain* This, char*& hash, int& size ) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialised atlas_Domain" );
    eckit::MD5 md5;
    This->hash( md5 );
    std::string s = md5.digest();
    size          = static_cast<int>( s.size() + 1 );
    hash          = new char[size];
    strcpy( hash, s.c_str() );
}
Domain::Spec* atlas__Domain__spec( const Domain* This ) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialised atlas_Domain" );
    return new Domain::Spec( This->spec() );
}

}  // extern "C"

}  // namespace domain
}  // namespace atlas
