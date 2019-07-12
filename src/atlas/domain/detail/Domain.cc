/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

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

}  // namespace domain
}  // namespace atlas
