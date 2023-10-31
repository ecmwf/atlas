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
    std::string type;
    config->get("type", type);
    auto redist = redistribution::detail::RedistributionImplFactory::build(type);
    redist->setup(fspace1, fspace2);
    return redist;
}

}


// ----------------------------------------------------------------------------

}  // namespace redistribution
}  // namespace atlas
