/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <cstring>

#include "RedistributeGeneric.h"
#include "RedistributionImpl.h"
#include "RedistributionInterface.h"

#include "atlas/redistribution/detail/RedistributionImplFactory.h"

namespace atlas {
namespace redistribution {

// ----------------------------------------------------------------------------
// Fortran interfaces
// ----------------------------------------------------------------------------

extern "C" {

detail::RedistributionImpl* atlas__Redistribution__new__config(
    const functionspace::FunctionSpaceImpl* fspace1, const functionspace::FunctionSpaceImpl* fspace2, 
    const eckit::Configuration* config) {
    ATLAS_ASSERT(config != nullptr);
    std::string type = detail::RedistributeGeneric::static_type();
    config->get("type", type);
    auto redist = redistribution::detail::RedistributionImplFactory::build(type);
    redist->setup(fspace1, fspace2);
    return redist;
}

void atlas__Redistribution__execute(
    const detail::RedistributionImpl* This, const Field* field_1, Field* field_2) {
    printf("redist execute");
    This->execute(*field_1, *field_2);
}

const FunctionSpace* atlas__Redistribution__source(
    const detail::RedistributionImpl* This) {
    return &This->source();
}

const FunctionSpace* atlas__Redistribution__target(
    const detail::RedistributionImpl* This) {
    return &This->target();
}

}


// ----------------------------------------------------------------------------

}  // namespace redistribution
}  // namespace atlas
