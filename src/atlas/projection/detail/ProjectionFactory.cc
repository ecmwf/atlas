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

#include "atlas/projection/Projection.h"
#include "atlas/projection/detail/ProjectionFactory.h"
#include "atlas/util/Config.h"

#include "atlas/projection/detail/CubedSphereEquiAnglProjection.h"
#include "atlas/projection/detail/CubedSphereEquiDistProjection.h"
#include "atlas/projection/detail/LambertAzimuthalEqualAreaProjection.h"
#include "atlas/projection/detail/LambertConformalConicProjection.h"
#include "atlas/projection/detail/LonLatProjection.h"
#include "atlas/projection/detail/MercatorProjection.h"
#include "atlas/projection/detail/SchmidtProjection.h"
#include "atlas/projection/detail/VariableResolutionProjection.h"

namespace atlas {
namespace projection {

//----------------------------------------------------------------------------------------------------------------------

void force_link() {
    static struct Link {
        Link() {
            ProjectionBuilder<detail::CubedSphereEquiAnglProjection>();
            ProjectionBuilder<detail::CubedSphereEquiDistProjection>();
            ProjectionBuilder<detail::LonLatProjection>();
            ProjectionBuilder<detail::RotatedLonLatProjection>();
            ProjectionBuilder<detail::SchmidtProjection>();
            ProjectionBuilder<detail::RotatedSchmidtProjection>();
            ProjectionBuilder<detail::MercatorProjection>();
            ProjectionBuilder<detail::RotatedMercatorProjection>();
            ProjectionBuilder<detail::VariableResolutionProjection>();
            ProjectionBuilder<detail::RotatedVariableResolutionProjection>();
            ProjectionBuilder<detail::LambertConformalConicProjection>();
            ProjectionBuilder<detail::LambertAzimuthalEqualAreaProjection>();
        }
    } link;
}

//----------------------------------------------------------------------------------------------------------------------

const Projection::Implementation* ProjectionFactory::build(const std::string& builder) {
    return build(builder, util::NoConfig());
}

const Projection::Implementation* ProjectionFactory::build(const std::string& builder,
                                                           const eckit::Parametrisation& param) {
    force_link();
    auto factory = get(builder);
    return factory->make(param);
}

//----------------------------------------------------------------------------------------------------------------------

}  // namespace projection
}  // namespace atlas
