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

#include "eckit/utils/Hash.h"

#include "atlas/domain/detail/DomainFactory.h"
#include "atlas/domain/detail/EmptyDomain.h"
#include "atlas/runtime/Exception.h"

namespace atlas {
namespace domain {

EmptyDomain::EmptyDomain() = default;

EmptyDomain::EmptyDomain( const eckit::Parametrisation& p ) {}

EmptyDomain::Spec EmptyDomain::spec() const {
    Spec domain_spec;
    domain_spec.set( "type", type() );
    return domain_spec;
}

void EmptyDomain::print( std::ostream& os ) const {
    os << "EmptyDomain";
}

void EmptyDomain::hash( eckit::Hash& h ) const {
    h.add( type() );
}

std::string EmptyDomain::units() const {
    ATLAS_NOTIMPLEMENTED;
}

namespace {
static DomainBuilder<EmptyDomain> register_builder( EmptyDomain::static_type() );
}

}  // namespace domain
}  // namespace atlas
