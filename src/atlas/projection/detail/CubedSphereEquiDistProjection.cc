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
  auto cubeNx = getCubeNx();
  auto tile1Lats = getLatArray();
  auto tile1Lons = getLonArray();

  // Cartesian coordinate starting points
  const double rsq3 = 1.0/sqrt(3.0);
  const double x0 = -rsq3;
  const double y0 = rsq3;
  const double z0 = -rsq3;

  // Equidistant
  const double dy = -2.0*rsq3/cubeNx;
  const double dz =  2.0*rsq3/cubeNx;

  double p1;
  double p2;
  double p3;
  double lonlat[2];

  for ( int ix = 0; ix < cubeNx+1; ix++ ) {
    for ( int iy = 0; iy < cubeNx+1; iy++ ) {
      // Grid points in cartesian coordinates
      double p1 = x0;
      double p2 = y0 + iy*dy;
      double p3 = z0 + ix*dz;

      ProjectionUtilities::cartesianToLatLon(p1, p2, p3, lonlat);

      tile1Lons(ix,iy) = lonlat[LON] - M_PI;
      tile1Lats(ix,iy) = lonlat[LAT];
    }
  }

  std::cout << "tile1Lons(0,0)"           << tile1Lons(0,0) << std::endl;
  std::cout << "tile1Lons(cubeNx,cubeNx)" << tile1Lons(cubeNx,cubeNx) << std::endl;
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
