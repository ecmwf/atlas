/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "CubedSphereEquiAnglProjection.h"

#include <cmath>

#include "eckit/config/Parametrisation.h"
#include "eckit/types/FloatCompare.h"
#include "eckit/utils/Hash.h"

#include "atlas/projection/detail/ProjectionFactory.h"
#include "atlas/runtime/Exception.h"
#include "atlas/util/Config.h"
#include "atlas/util/Constants.h"
#include "atlas/util/CoordinateEnums.h"
#include "atlas/util/Earth.h"

namespace atlas {
namespace projection {
namespace detail {

CubedSphereEquiAnglProjection::CubedSphereEquiAnglProjection( const eckit::Parametrisation& params ) {}

void CubedSphereEquiAnglProjection::lonlat2xy( double crd[] ) const {}

void CubedSphereEquiAnglProjection::xy2lonlat( double crd[] ) const {}

CubedSphereEquiAnglProjection::Spec CubedSphereEquiAnglProjection::spec() const {}

void CubedSphereEquiAnglProjection::hash( eckit::Hash& h ) const {}

namespace {
static ProjectionBuilder<CubedSphereEquiAnglProjection>
       register_1( CubedSphereEquiAnglProjection::static_type() );
}

}  // namespace detail
}  // namespace projection
}  // namespace atlas
