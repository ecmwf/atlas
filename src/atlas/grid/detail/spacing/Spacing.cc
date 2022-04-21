/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "Spacing.h"

#include "eckit/config/Parametrisation.h"

#include "atlas/grid/detail/spacing/SpacingFactory.h"
#include "atlas/runtime/Exception.h"

namespace atlas {
namespace grid {
namespace spacing {

const Spacing* Spacing::create(const eckit::Parametrisation& params) {
    std::string spacingType;
    if (not params.get("type", spacingType)) {
        throw_Exception("type missing in configuration", Here());
    }
    return SpacingFactory::build(spacingType, params);
}

}  // namespace spacing
}  // namespace grid
}  // namespace atlas
