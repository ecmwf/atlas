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

#include "atlas/domain/detail/Domain.h"
#include "atlas/domain/detail/DomainFactory.h"
#include "atlas/domain/detail/EmptyDomain.h"
#include "atlas/domain/detail/GlobalDomain.h"
#include "atlas/domain/detail/RectangularDomain.h"
#include "atlas/domain/detail/ZonalBandDomain.h"

namespace atlas {
namespace domain {

//----------------------------------------------------------------------------------------------------------------------

namespace {
void force_link() {
    static struct Link {
        Link() {
            DomainBuilder<GlobalDomain>();
            DomainBuilder<EmptyDomain>();
            DomainBuilder<RectangularDomain>();
            DomainBuilder<ZonalBandDomain>();
        }
    } link;
}
}  // namespace

//----------------------------------------------------------------------------------------------------------------------

const Domain* DomainFactory::build(const std::string& builder) {
    return build(builder, util::NoConfig());
}

const Domain* DomainFactory::build(const std::string& builder, const eckit::Parametrisation& param) {
    force_link();
    auto factory = get(builder);
    return factory->make(param);
}

//----------------------------------------------------------------------------------------------------------------------

}  // namespace domain
}  // namespace atlas
