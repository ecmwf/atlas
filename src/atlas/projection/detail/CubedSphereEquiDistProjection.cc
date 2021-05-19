/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "CubedSphereEquiDistProjection.h"

#include <cmath>
#include <iostream>

#include "eckit/config/Parametrisation.h"
#include "eckit/types/FloatCompare.h"
#include "eckit/utils/Hash.h"

#include "atlas/projection/detail/ProjectionFactory.h"
#include "atlas/projection/detail/ProjectionUtilities.h"
#include "atlas/runtime/Exception.h"
#include "atlas/util/Config.h"
#include "atlas/util/Constants.h"
#include "atlas/util/CoordinateEnums.h"
#include "atlas/util/Earth.h"

namespace atlas {
namespace projection {
namespace detail {

using util::Constants;

// -------------------------------------------------------------------------------------------------

CubedSphereEquiDistProjection::CubedSphereEquiDistProjection( const eckit::Parametrisation& params )
    : CubedSphereProjectionBase(params) {
}

// -------------------------------------------------------------------------------------------------

void CubedSphereEquiDistProjection::lonlat2xy( double crd[] ) const {

    std::cout << "lonlat2xy start : lonlat = " << crd[0] << " " << crd[1] << std::endl;

    double xyz[3];
    double ab[2]; // alpha-beta coordinate

    // convert degrees to radians
    crd[0] *= M_PI/180.;
    crd[1] *= M_PI/180.;

    // To [-pi/4, 7/8 *pi)
   if (crd[LON] >= 1.75 * M_PI) crd[LON] += -2.*M_PI;

    // find tile which this lonlat is linked to
    idx_t t =  CubedSphereProjectionBase::tileFromLonLat(crd);

    ProjectionUtilities::sphericalToCartesian(crd, xyz,false, true);
    tileRotateInverse.at(t)(xyz);

    //now should be tile 0 - now calculate (alpha, beta) in radians.
    // should be between - pi/4 and pi/4
    ab[0] =   M_PI_4 * xyz[YY] / xyz[XX] ;
    ab[1] = - M_PI_4 * xyz[ZZ] / xyz[XX];

    std::cout << "lonlat2xy xyz ab : "
       << xyz[0] << " " << xyz[1]  << " " << xyz[2] << " "
       << ab[0] << " " << ab[1] << std::endl;

    CubedSphereProjectionBase::alphabetatt2xy(t, ab, crd);

    std::cout << "lonlat2xy end : xy = " << crd[0] << " " << crd[1] << std::endl;

}


// -------------------------------------------------------------------------------------------------

void CubedSphereEquiDistProjection::xy2lonlat( double crd[] ) const {

    const double rsq3 = 1.0/sqrt(3.0);
    double xyz[3];
    double ab[2]; // alpha-beta coordinate
    idx_t t;  // tile index

    // calculate xy (in degrees) to alpha beta (in radians) and t - tile index.
    CubedSphereProjectionBase::xy2alphabetat(crd, t, ab);

    std::cout << "xy2lonlat:: crd t ab  : "  << crd[0] << " " << crd[1] << " " << t << " " << ab[0] << " " << ab[1] << std::endl;

    xyz[0] = -rsq3;
    xyz[1] = -rsq3 * ab[0] / M_PI_4;
    xyz[2] = -rsq3 * ab[1] / M_PI_4;

    ProjectionUtilities::cartesianToSpherical(xyz, crd, false);

    if (crd[LON] < 0.0) crd[LON] += 2.0*M_PI;
    crd[LON] = crd[LON] - M_PI;

    std::cout << "xy2lonlat:: lonlat before rotation : "  << crd[0] << " " << crd[1]  << std::endl;

    // Convert to cartesian
    ProjectionUtilities::sphericalToCartesian(crd, xyz, false, true);

    // Perform tile specific rotation
    tileRotate.at(t)(xyz);

    // Back to latlon
    ProjectionUtilities::cartesianToSpherical(xyz, crd, false);

    // Shift longitude
    /*
    if (shiftLon_ != 0.0) {
      crd[LON] = crd[LON] + shiftLon_*atlas::util::Constants::degreesToRadians();
      if (crd[LON] < -M_PI) {crd[LON] =  2*M_PI + crd[LON];}
      if (crd[LON] >  M_PI) {crd[LON] = -2*M_PI + crd[LON];}
    }
    */
    // To 0, 360
    if (crd[LON] < 0.0) crd[LON] = 2*M_PI + crd[LON];

    // longitude does not make sense at the poles - set to 0.
    if ( std::abs(std::abs(crd[LAT]) - M_PI_2) < 1e-13) crd[LON] = 0.;

    crd[LON] *= Constants::radiansToDegrees();
    crd[LAT] *= Constants::radiansToDegrees();

    std::cout << "end of equidistant xy2lonlat lonlat = " <<  crd[LON] << " " << crd[LAT] << std::endl;
}


// -------------------------------------------------------------------------------------------------

ProjectionImpl::Jacobian CubedSphereEquiDistProjection::jacobian(const PointLonLat& ) const {
    ATLAS_NOTIMPLEMENTED;
}

// -------------------------------------------------------------------------------------------------

// -------------------------------------------------------------------------------------------------

CubedSphereEquiDistProjection::Spec CubedSphereEquiDistProjection::spec() const {
    // Fill projection specification
    Spec proj;
    proj.set( "type", static_type() );
    return proj;
}

// -------------------------------------------------------------------------------------------------

void CubedSphereEquiDistProjection::hash( eckit::Hash& h ) const {
    // Add to hash
    h.add( static_type() );
    CubedSphereProjectionBase::hash(h);
}

// -------------------------------------------------------------------------------------------------

namespace {
static ProjectionBuilder<CubedSphereEquiDistProjection>
register_1( CubedSphereEquiDistProjection::static_type() );
}

}  // namespace detail
}  // namespace projection
}  // namespace atlas
