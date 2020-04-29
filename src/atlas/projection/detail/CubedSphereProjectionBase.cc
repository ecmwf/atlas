/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <cmath>

#include "CubedSphereProjectionBase.h"

#include "atlas/projection/detail/ProjectionUtilities.h"
#include "atlas/util/Constants.h"
#include "atlas/util/CoordinateEnums.h"

namespace atlas {
namespace projection {
namespace detail {

// -------------------------------------------------------------------------------------------------

CubedSphereProjectionBase::CubedSphereProjectionBase( const eckit::Parametrisation& params )
                                                            : tile1LonsArray_(), tile1LatsArray_() {
  ATLAS_TRACE( "CubedSphereProjectionBase::CubedSphereProjectionBase" );
  // Get cube sphere face dimension
  params.get("CubeNx", cubeNx_);

  // Arrays to hold projection for (0,0) centered tile
  tile1LonsArray_.reset(new ArrayLatLon_(cubeNx_+1, cubeNx_+1));
  tile1LatsArray_.reset(new ArrayLatLon_(cubeNx_+1, cubeNx_+1));
}

// -------------------------------------------------------------------------------------------------

void CubedSphereProjectionBase::lonlat2xy( double llxytl[] ) const {

  double lonlat[2];
  double lonlat_in[2];
  double xyz[3];

  // Get full tile 1 lonlat arrays
  auto tile1Lats = getLatArray();
  auto tile1Lons = getLonArray();

  // Get lonlat point array
  lonlat_in[LON] = llxytl[LON];
  lonlat_in[LAT] = llxytl[LAT];

  // To -180, 180
  if (lonlat_in[LON] > 180.0) {
    lonlat_in[LON] = lonlat_in[LON] - 360.0;
  }

  // To radians
  lonlat_in[LON] = lonlat_in[LON] * atlas::util::Constants::degreesToRadians();
  lonlat_in[LAT] = lonlat_in[LAT] * atlas::util::Constants::degreesToRadians();

  // Loop over tiles
  bool found = false;
  for (int t = 0; t < 6; t++) {

    if (!found) {

      // Copy before overwriting
      std::copy(lonlat_in, lonlat_in+2, lonlat);

      // Convert to cartesian
      ProjectionUtilities::sphericalToCartesian(lonlat, xyz);

      // Perform tile specific rotations
      tileRotateInverse.at(t)(xyz);

      // Back to latlon
      ProjectionUtilities::cartesianToSpherical(xyz, lonlat);

      // Loop over tile 1 array
      for (int ix = 0; ix < cubeNx_; ix++) {
        for (int iy = 0; iy < cubeNx_; iy++) {
          double absLonDiff = std::abs(tile1Lons(ix, iy) - lonlat[LON]);
          double absLatDiff = std::abs(tile1Lats(ix, iy) - lonlat[LAT]);
          // Check if this rotation landed the point on tile 1
          if (!found && absLonDiff < 1.0e-15 && absLatDiff < 1.0e-15) {
            llxytl[2+XX] = ix;
            llxytl[2+YY] = iy;
            llxytl[4] = t;
            found = true;
          }
        }
      }
    }
  }
  ATLAS_ASSERT(found, "CubedSphereProjectionBase::lonlat2xy, LonLat not found");
}

// -------------------------------------------------------------------------------------------------

void CubedSphereProjectionBase::xy2lonlat( double xytll[] ) const {

  double lonlat[2];
  double xyz[3];
  double angle;

  // Get lat/lon for this index but on tile 1
  auto tile1Lats = getLatArray();
  auto tile1Lons = getLonArray();
  lonlat[LON] = tile1Lons(static_cast<int>(xytll[XX]), static_cast<int>(xytll[YY]));
  lonlat[LAT] = tile1Lats(static_cast<int>(xytll[XX]), static_cast<int>(xytll[YY]));

  // Convert to cartesian
  ProjectionUtilities::sphericalToCartesian(lonlat, xyz);

  // Perform tile specific rotation
  tileRotate.at(xytll[2])(xyz);

  // Back to latlon
  ProjectionUtilities::cartesianToSpherical(xyz, lonlat);

  // Fill outputs and covert to degrees
  xytll[3+LON] = lonlat[LON] * atlas::util::Constants::radiansToDegrees();
  xytll[3+LAT] = lonlat[LAT] * atlas::util::Constants::radiansToDegrees();

  // To 0, 360
  if (xytll[3+LON] < 0.0) {
    xytll[3+LON] = 360.0 + xytll[3+LON];
  }

}

// -------------------------------------------------------------------------------------------------

void CubedSphereProjectionBase::tile1Rotate( double xyz[] ) const {
  //  Face 1, no rotation.
}

// -------------------------------------------------------------------------------------------------

void CubedSphereProjectionBase::tile2Rotate( double xyz[] ) const {
  //  Face 2: rotate -90.0 degrees about z axis
  double angle;
  angle = -M_PI / 2.0;
  ProjectionUtilities::rotate3dZ(angle, xyz);
}

// -------------------------------------------------------------------------------------------------

void CubedSphereProjectionBase::tile3Rotate( double xyz[] ) const {
  //  Face 3: rotate -90.0 degrees about z axis
  //          rotate  90.0 degrees about x axis
  double angle;
  angle = -M_PI / 2.0;
  ProjectionUtilities::rotate3dZ(angle, xyz);
  angle = M_PI / 2.0;
  ProjectionUtilities::rotate3dX(angle, xyz);
}

// -------------------------------------------------------------------------------------------------

void CubedSphereProjectionBase::tile4Rotate( double xyz[] ) const {
  //  Face 4: rotate -180.0 degrees about z axis
  //          rotate   90.0 degrees about x axis
  double angle;
  angle = -M_PI;
  ProjectionUtilities::rotate3dZ(angle, xyz);
  angle = M_PI / 2.0;
  ProjectionUtilities::rotate3dX(angle, xyz);
}

// -------------------------------------------------------------------------------------------------

void CubedSphereProjectionBase::tile5Rotate( double xyz[] ) const {
  //  Face 5: rotate 90.0 degrees about z axis
  //          rotate 90.0 degrees about y axis
  double angle;
  angle = M_PI / 2.0;
  ProjectionUtilities::rotate3dZ(angle, xyz);
  angle = M_PI / 2.0;
  ProjectionUtilities::rotate3dY(angle, xyz);
}

// -------------------------------------------------------------------------------------------------

void CubedSphereProjectionBase::tile6Rotate( double xyz[] ) const {
  //  Face 6: rotate 90.0 degrees about y axis
  double angle;
  angle = M_PI / 2.0;
  ProjectionUtilities::rotate3dY(angle, xyz);
}

// -------------------------------------------------------------------------------------------------

void CubedSphereProjectionBase::tile1RotateInverse( double llxyt[] ) const {
  //  Face 1, no rotation.
}

// -------------------------------------------------------------------------------------------------

void CubedSphereProjectionBase::tile2RotateInverse( double xyz[] ) const {
  //  Face 2: rotate 90.0 degrees about z axis
  double angle;
  angle = M_PI / 2.0;
  ProjectionUtilities::rotate3dZ(angle, xyz);
}

// -------------------------------------------------------------------------------------------------

void CubedSphereProjectionBase::tile3RotateInverse( double xyz[] ) const {
  //  Face 3: rotate  90.0 degrees about x axis
  //          rotate -90.0 degrees about z axis
  double angle;
  angle = -M_PI / 2.0;
  ProjectionUtilities::rotate3dX(angle, xyz);
  angle = M_PI / 2.0;
  ProjectionUtilities::rotate3dZ(angle, xyz);
}

// -------------------------------------------------------------------------------------------------

void CubedSphereProjectionBase::tile4RotateInverse( double xyz[] ) const {
  //  Face 4: rotate   90.0 degrees about x axis
  //          rotate -180.0 degrees about z axis
  double angle;
  angle = -M_PI / 2.0;
  ProjectionUtilities::rotate3dX(angle, xyz);
  angle = M_PI;
  ProjectionUtilities::rotate3dZ(angle, xyz);
}

// -------------------------------------------------------------------------------------------------

void CubedSphereProjectionBase::tile5RotateInverse( double xyz[] ) const {
  //  Face 5: rotate -90.0 degrees about y axis
  //          rotate -90.0 degrees about z axis
  double angle;
  angle = -M_PI / 2.0;
  ProjectionUtilities::rotate3dY(angle, xyz);
  angle = -M_PI / 2.0;
  ProjectionUtilities::rotate3dZ(angle, xyz);
}

// -------------------------------------------------------------------------------------------------

void CubedSphereProjectionBase::tile6RotateInverse( double xyz[] ) const {
  //  Face 6: rotate -90.0 degrees about y axis
  double angle;
  angle = -M_PI / 2.0;
  ProjectionUtilities::rotate3dY(angle, xyz);
}

// -------------------------------------------------------------------------------------------------

}  // namespace detail
}  // namespace projection
}  // namespace atlas
