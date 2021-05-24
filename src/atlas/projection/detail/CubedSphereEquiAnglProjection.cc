/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "CubedSphereEquiAnglProjection.h"

#include <cmath>
#include <iostream>

#include "eckit/config/Parametrisation.h"
#include "eckit/types/FloatCompare.h"
#include "eckit/utils/Hash.h"

#include "atlas/projection/detail/ProjectionFactory.h"
#include "atlas/projection/detail/ProjectionUtilities.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"
#include "atlas/util/Config.h"
#include "atlas/util/Constants.h"
#include "atlas/util/CoordinateEnums.h"

namespace atlas {
namespace projection {
namespace detail {

// -------------------------------------------------------------------------------------------------

CubedSphereEquiAnglProjection::CubedSphereEquiAnglProjection( const eckit::Parametrisation& params )
    : CubedSphereProjectionBase(params) {

}

// -------------------------------------------------------------------------------------------------

void CubedSphereEquiAnglProjection::lonlat2xy( double crd[] ) const {

    Log::info() << "lonlat2xy start : lonlat = " << crd[LON] << " " << crd[LAT] << std::endl;

    idx_t t;
    double ab[2]; // alpha-beta coordinate
    double xyz[3]; // on Cartesian grid

    CubedSphereProjectionBase::lonlat2xypre(crd, t, xyz);

    // should be between - pi/4 and pi/4
    // now calculate (alpha, beta) in radians.
    ab[0] = std::atan(xyz[YY]/xyz[XX]);
    ab[1] = std::atan(-xyz[ZZ]/xyz[XX]);  // I think the minus is here due to the
    // left coordinate system

    Log::debug() << "lonlat2xy xyz ab : "
      << xyz[XX] << " " << xyz[YY] << " " << xyz[ZZ] << " "
      << ab[LON] << " " << ab[LAT] << std::endl;

    CubedSphereProjectionBase::alphabetatt2xy(t, ab, crd);

    Log::debug() << "lonlat2xy end : xy = " << crd[LON] << " " << crd[LAT] << std::endl;

}

// -------------------------------------------------------------------------------------------------
// input should be Willems xy coordinate in degrees
//
void CubedSphereEquiAnglProjection::xy2lonlat( double crd[] ) const {

    Log::info() << "xy2lonlat start xy = " << crd[LON] << " " << crd[LAT] <<std::endl;

    const double rsq3 = 1.0/std::sqrt(3.0);
    double xyz[3];
    double ab[2]; // alpha-beta coordinate
    idx_t t;  // tile index

    // calculate xy (in degrees) to alpha beta (in radians) and t - tile index.
    CubedSphereProjectionBase::xy2alphabetat(crd, t, ab);

    Log::debug() << "xy2lonlat:: crd t ab  : "  << crd[0] << " " << crd[1] << " " << t << " " << ab[0] << " " << ab[1] << std::endl;

    xyz[0] = -rsq3;
    xyz[1] = -rsq3*tan(ab[0]);
    xyz[2] = -rsq3*tan(ab[1]);

    CubedSphereProjectionBase::xy2lonlatpost(xyz, t, crd);

    Log::info() << "end of equiangular xy2lonlat lonlat = " <<  crd[LON] << " " << crd[LAT] << std::endl;
}

// -------------------------------------------------------------------------------------------------

ProjectionImpl::Jacobian CubedSphereEquiAnglProjection::jacobian(const PointLonLat& ) const {
    ATLAS_NOTIMPLEMENTED;
}

// -------------------------------------------------------------------------------------------------

CubedSphereEquiAnglProjection::Spec CubedSphereEquiAnglProjection::spec() const {
    // Fill projection specification
    Spec proj;
    proj.set( "type", static_type() );
    return proj;
}

// -------------------------------------------------------------------------------------------------

void CubedSphereEquiAnglProjection::hash( eckit::Hash& h ) const {
    // Add to hash
    h.add( static_type() );
    CubedSphereProjectionBase::hash(h);
}

// -------------------------------------------------------------------------------------------------

namespace {
static ProjectionBuilder<CubedSphereEquiAnglProjection>
register_1( CubedSphereEquiAnglProjection::static_type() );
}

}  // namespace detail
}  // namespace projection
}  // namespace atlas
