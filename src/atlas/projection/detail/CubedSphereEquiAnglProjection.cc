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
#include "atlas/util/Config.h"
#include "atlas/util/Constants.h"
#include "atlas/util/CoordinateEnums.h"
#include "atlas/util/Earth.h"


namespace atlas {
namespace projection {
namespace detail {

// -------------------------------------------------------------------------------------------------

CubedSphereEquiAnglProjection::CubedSphereEquiAnglProjection( const eckit::Parametrisation& params )
                                                               : CubedSphereProjectionBase(params) {


 /*
  // Get Base data
  const auto cubeNx = getCubeNx();
  auto tile1Lats = getLatArray();
  auto tile1Lons = getLonArray();

  // Cartesian coordinate starting points
  const double rsq3 = 1.0/sqrt(3.0);

  // Equiangular
  const double dp = 0.5*M_PI/static_cast<double>(cubeNx);

  double xyz[3];
  double lonlat[2];

  for ( int ix = 0; ix < cubeNx+1; ix++ ) {
    for ( int iy = 0; iy < cubeNx+1; iy++ ) {
      // Grid points in cartesian coordinates
      xyz[XX] = -rsq3;
      xyz[YY] = -rsq3*tan(-0.25*M_PI+static_cast<double>(ix)*dp);
      xyz[ZZ] = -rsq3*tan(-0.25*M_PI+static_cast<double>(iy)*dp);

      ProjectionUtilities::cartesianToSpherical(xyz, lonlat);

      if (lonlat[LON] < 0.0) {
        lonlat[LON] += 2.0*M_PI;
      }

      tile1Lons(ix, iy) = lonlat[LON] - M_PI;
      tile1Lats(ix, iy) = lonlat[LAT];
    }
  }
  */
}

// -------------------------------------------------------------------------------------------------

void CubedSphereEquiAnglProjection::lonlat2xy( double crd[] ) const {

   std::cout << "lonlat2xy start : lonlat = " << crd[0] << " " << crd[1] << std::endl;

   if (std::abs(crd[LON]) < 1e-14) crd[LON] = 0.;
   if (std::abs(crd[LAT]) < 1e-14) crd[LAT] = 0.;

   double xyz[3];
   double ab[2]; // alpha-beta coordinate

   // convert degrees to radians
   crd[0] *= M_PI/180.;
   crd[1] *= M_PI/180.;

   // To [-pi/4, 7/8 *pi]
  if (crd[LON] >= 1.75 * M_PI) {
    crd[LON] += -2.*M_PI;
  }

   // find tile which this lonlat is linked to
   // works [0,2pi]
   idx_t t = CubedSphereProjectionBase::tileFromLonLat(crd);

   /*
   if (crd[LON] < 0.0) {
     crd[LON] += 2.0*M_PI;
   }
   crd[LON] = crd[LON] - M_PI;
   */

   /*

   */

   ProjectionUtilities::sphericalToCartesian(crd, xyz, false, true);
   tileRotateInverse.at(t)(xyz);

   //now should be tile 0 - now calculate (alpha, beta) in radians.
   // should be between - pi/4 and pi/4
   ab[0] = atan(xyz[1]/xyz[0]);
   ab[1] = atan(-xyz[2]/xyz[0]);  // I think the minus is here due to the
                                  // left coordinate system

   std::cout << "lonlat2xy xyz ab : "
      << xyz[0] << " " << xyz[1]  << " " << xyz[2] << " "
      << ab[0] << " " << ab[1] << std::endl;

   CubedSphereProjectionBase::alphabetatt2xy(t, ab, crd);

   std::cout << "lonlat2xy end : xy = " << crd[0] << " " << crd[1] << std::endl;


}

// -------------------------------------------------------------------------------------------------
// input should be Willems xy coordinate in degrees
//
void CubedSphereEquiAnglProjection::xy2lonlat( double crd[] ) const {

    std::cout << "xy2lonlat start xy = " << crd[LON] << " " << crd[LAT] <<std::endl;

    const double rsq3 = 1.0/sqrt(3.0);
    double xyz[3];
    double ab[2]; // alpha-beta coordinate
    idx_t t;  // tile index
    double lonlat[2];

    // calculate xy (in degrees) to alpha beta (in radians) and t - tile index.
    CubedSphereProjectionBase::xy2alphabetat(crd, t, ab);

    std::cout << "xy2lonlat:: crd t ab  : "  << crd[0] << " " << crd[1] << " " << t << " " << ab[0] << " " << ab[1] << std::endl;

    xyz[0] = -rsq3;
    xyz[1] = -rsq3*tan(ab[0]);
    xyz[2] = -rsq3*tan(ab[1]);

    ProjectionUtilities::cartesianToSpherical(xyz, lonlat, false);

    if (lonlat[LON] < 0.0) {
      lonlat[LON] += 2.0*M_PI;
    }
    lonlat[LON] = lonlat[LON] - M_PI;

    std::cout << "xy2lonlat:: lonlat before rotation : "  << lonlat[0] << " " << lonlat[1]  << std::endl;

    // Convert to cartesian
    ProjectionUtilities::sphericalToCartesian(lonlat, xyz, false, true);

    // Perform tile specific rotation
    tileRotate.at(t)(xyz);

    // Back to latlon
    ProjectionUtilities::cartesianToSpherical(xyz, lonlat, false);

    // Shift longitude
    /*
    if (shiftLon_ != 0.0) {
      lonlat[LON] = lonlat[LON] + shiftLon_*atlas::util::Constants::degreesToRadians();
      if (lonlat[LON] < -M_PI) {lonlat[LON] =  2*M_PI + lonlat[LON];}
      if (lonlat[LON] >  M_PI) {lonlat[LON] = -2*M_PI + lonlat[LON];}
    }
    */
    // To 0, 360
    if (lonlat[LON] < 0.0) {
      lonlat[LON] = 2.*M_PI + lonlat[LON];
    }

    // longitude does not make sense at the poles - set to 0.
    if ( std::abs(std::abs(lonlat[LAT]) - M_PI_2) < 1e-13) lonlat[LON] = 0.;


    crd[LON] = lonlat[LON] * 180.0 / M_PI;
    crd[LAT] = lonlat[LAT] * 180.0 / M_PI;


    std::cout << "end of xy2lonlat lonlat = " <<  crd[LON] << " " << crd[LAT] << std::endl;
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
