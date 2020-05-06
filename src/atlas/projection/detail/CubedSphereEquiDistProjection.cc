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

// -------------------------------------------------------------------------------------------------

CubedSphereEquiDistProjection::CubedSphereEquiDistProjection( const eckit::Parametrisation& params )
                                                               : CubedSphereProjectionBase(params) {
  // Get Base data
  const auto cubeNx = getCubeNx();
  auto tile1Lats = getLatArray();
  auto tile1Lons = getLonArray();

  // Cartesian coordinate starting points
  const double rsq3 = 1.0/sqrt(3.0);

  // Equidistant
  const double dy = -2.0*rsq3/static_cast<double>(cubeNx);
  const double dz =  2.0*rsq3/static_cast<double>(cubeNx);

  double xyz[3];
  double lonlat[2];

  for ( int ix = 0; ix < cubeNx+1; ix++ ) {
    for ( int iy = 0; iy < cubeNx+1; iy++ ) {
      // Grid points in cartesian coordinates
      xyz[XX] = -rsq3;
      xyz[YY] =  rsq3 + static_cast<double>(ix)*dy;
      xyz[ZZ] =  rsq3 - static_cast<double>(iy)*dz;

      ProjectionUtilities::cartesianToSpherical(xyz, lonlat);

      if (lonlat[LON] < 0.0) {
        lonlat[LON] += 2.0*M_PI;
      }

      tile1Lons(ix, iy) = lonlat[LON] - M_PI;
      tile1Lats(ix, iy) = lonlat[LAT];
    }
  }
}

// -------------------------------------------------------------------------------------------------

void CubedSphereEquiDistProjection::lonlat2xy( double crd[] ) const {
  CubedSphereProjectionBase::lonlat2xy(crd);
}

// -------------------------------------------------------------------------------------------------

void CubedSphereEquiDistProjection::xy2lonlat( double crd[] ) const {
  CubedSphereProjectionBase::xy2lonlat(crd);
}

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
}

// -------------------------------------------------------------------------------------------------

namespace {
static ProjectionBuilder<CubedSphereEquiDistProjection>
       register_1( CubedSphereEquiDistProjection::static_type() );
}

}  // namespace detail
}  // namespace projection
}  // namespace atlas
