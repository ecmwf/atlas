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

#include "atlas/functionspace/FunctionSpace.h"
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
    FunctionSpace fs1(fspace1);
    FunctionSpace fs2(fspace2);
    redist->setup(fs1, fs2);
    return redist;
}

void atlas__Redistribution__execute(
    const detail::RedistributionImpl* This, const field::FieldImpl* field_1, field::FieldImpl* field_2) {
    Field f1(field_1);
    Field f2(field_2);
    This->execute(f1, f2);
}

const functionspace::FunctionSpaceImpl* atlas__Redistribution__source(
    const detail::RedistributionImpl* This) {
    return This->source().get();
}

const functionspace::FunctionSpaceImpl* atlas__Redistribution__target(
    const detail::RedistributionImpl* This) {
    return This->target().get();
}

}


// ----------------------------------------------------------------------------

}  // namespace redistribution
}  // namespace atlas
