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

  // Shift projection by a longitude
  shiftLon_ = 0.0;
  if (params.has("ShiftLon")) {
    params.get("ShiftLon", shiftLon_);
    ATLAS_ASSERT(shiftLon_ <= 90.0,  "ShiftLon should be <= 90.0 degrees");
    ATLAS_ASSERT(shiftLon_ >= -90.0, "ShiftLon should be >= -90.0 degrees");
  }

  // Apply a Schmidt transform
  doSchmidt_ = false;
  stretchFac_ = 0.0;
  targetLon_ = 0.0;
  targetLat_ = 0.0;
  if (params.has("DoSchmidt")) {
    params.get("DoSchmidt", doSchmidt_);
    if (doSchmidt_) {
      params.get("StretchFac", stretchFac_);
      params.get("TargetLon", targetLon_);
      params.get("TargetLat", targetLat_);
    }
  }

  // Arrays to hold projection for (0,0) centered tile
  tile1LonsArray_.reset(new ArrayLatLon_(cubeNx_+1, cubeNx_+1));
  tile1LatsArray_.reset(new ArrayLatLon_(cubeNx_+1, cubeNx_+1));
}

// -------------------------------------------------------------------------------------------------

void CubedSphereProjectionBase::hash( eckit::Hash& h ) const {
  // Add stretching options to hash
  h.add(shiftLon_);
  h.add(doSchmidt_);
  h.add(stretchFac_);
  h.add(targetLon_);
  h.add(targetLat_);
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
      tileRotateInverse.at(std::size_t(t))(xyz);

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

void CubedSphereProjectionBase::xy2lonlat( double crd[] ) const {

    // For now I think I will write the equiangular version here first:
    // 1) get tile.
    // 2) xy2alphabeta
    // 3) alphabeta2lonlat.
    //  look at Ronchi first.


  // double lonlat[2];
 // double xyz[3];

  // Get lat/lon for this index but on tile 1
  // auto tile1Lats = getLatArray();
  // auto tile1Lons = getLonArray();
  // lonlat[LON] = tile1Lons(static_cast<int>(xytll[XX]), static_cast<int>(xytll[YY]));
  // lonlat[LAT] = tile1Lats(static_cast<int>(xytll[XX]), static_cast<int>(xytll[YY]));

  //cmw - we require a transformation here that converts Willem's xy to lonlat.
  // the former method stored the lon lats for tile 1 -> converted to cartesian -> rotate grid
  //         -> converted to spherical.

  // xy2alphabetat
  // set xyz


   /*
   const double rsq3 = 1.0/sqrt(3.0);
   double xyz[3];
   xyz[XX] = -rsq3;
   xyz[YY] = -rsq3*tan(-0.25*M_PI+static_cast<double>(ix)*dp);
   xyz[ZZ] = -rsq3*tan(-0.25*M_PI+static_cast<double>(iy)*dp);

   ProjectionUtilities::cartesianToSpherical(xyz, lonlat);

   if (lonlat[LON] < 0.0) {
      lonlat[LON] += 2.0*M_PI;
   }

    tile1Lons(ix, iy) = lonlat[LON] - M_PI;
    tile1Lats(ix, iy) = lonlat[LAT];
  */

  /*

  // Convert to cartesian
  ProjectionUtilities::sphericalToCartesian(lonlat, xyz);

  // Perform tile specific rotation
  tileRotate.at(xytll[2])(xyz);

  // Back to latlon
  ProjectionUtilities::cartesianToSpherical(xyz, lonlat);

  // Shift longitude
  if (shiftLon_ != 0.0) {
    lonlat[LON] = lonlat[LON] + shiftLon_*atlas::util::Constants::degreesToRadians();
    if (lonlat[LON] < -M_PI) {lonlat[LON] =  2*M_PI + lonlat[LON];}
    if (lonlat[LON] >  M_PI) {lonlat[LON] = -2*M_PI + lonlat[LON];}
  }

  // To 0, 360
  if (lonlat[LON] < 0.0) {
    lonlat[LON] = 2*M_PI + lonlat[LON];
  }

  // Schmidt transform
  if (doSchmidt_) {
    this->schmidtTransform( stretchFac_, targetLon_*atlas::util::Constants::degreesToRadians(),
                            targetLat_*atlas::util::Constants::degreesToRadians(), lonlat);
  }

  // Fill outputs and covert to degrees
  xytll[3+LON] = lonlat[LON] * atlas::util::Constants::radiansToDegrees();
  xytll[3+LAT] = lonlat[LAT] * atlas::util::Constants::radiansToDegrees();

  */

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
  std::cout << "tile4 Rotate after 3dZ" << xyz[0]  << " " << xyz[1] << " " << xyz[2] << std::endl;
 // angle = M_PI / 2.0;
//  ProjectionUtilities::rotate3dX(angle, xyz);
  std::cout << "tile4 Rotate after 3dX" << xyz[0]  << " " << xyz[1] << " " << xyz[2] << std::endl;
}

// -------------------------------------------------------------------------------------------------

void CubedSphereProjectionBase::tile5Rotate( double xyz[] ) const {
  //  Face 5: rotate 90.0 degrees about z axis
  //          rotate 90.0 degrees about y axis
  double angle;
  angle = M_PI / 2.0;
  ProjectionUtilities::rotate3dZ(angle, xyz);
 // angle = M_PI / 2.0;
 // ProjectionUtilities::rotate3dY(angle, xyz);
}

// -------------------------------------------------------------------------------------------------

void CubedSphereProjectionBase::tile6Rotate( double xyz[] ) const {
  //  Face 6: rotate 90.0 degrees about y axis
  // double angle;
  // angle = M_PI / 2.0;
  // ProjectionUtilities::rotate3dY(angle, xyz);

  double angle;
// angle = - 3. * M_PI / 2.;
//  ProjectionUtilities::rotate3dZ(angle, xyz);
  angle =  M_PI /2. ;
  ProjectionUtilities::rotate3dY(angle, xyz);

  angle =  M_PI / 2.;
  ProjectionUtilities::rotate3dZ(angle, xyz);
}

// -------------------------------------------------------------------------------------------------

void CubedSphereProjectionBase::tile1RotateInverse( double xyt[] ) const {
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
//  angle = -M_PI / 2.0;
//  ProjectionUtilities::rotate3dX(angle, xyz);
  angle = M_PI;
  ProjectionUtilities::rotate3dZ(angle, xyz);
}

// -------------------------------------------------------------------------------------------------

void CubedSphereProjectionBase::tile5RotateInverse( double xyz[] ) const {
  //  Face 5: rotate -90.0 degrees about y axis
  //          rotate -90.0 degrees about z axis
  double angle;
//  angle = -M_PI / 2.0;
//  ProjectionUtilities::rotate3dY(angle, xyz);
  angle = - M_PI / 2.0;
  ProjectionUtilities::rotate3dZ(angle, xyz);
}

// -------------------------------------------------------------------------------------------------

void CubedSphereProjectionBase::tile6RotateInverse( double xyz[] ) const {
  //  Face 6: rotate -90.0 degrees about y axis
  double angle;
  angle = -M_PI / 2.;
  ProjectionUtilities::rotate3dZ(angle, xyz);
  angle = -M_PI / 2.;
  ProjectionUtilities::rotate3dY(angle, xyz);

}

// -------------------------------------------------------------------------------------------------

void CubedSphereProjectionBase::schmidtTransform( double stretchFac, double targetLon,
                                                  double targetLat, double lonlat[]) const {


  double c2p1 = 1.0 + stretchFac*stretchFac;
  double c2m1 = 1.0 - stretchFac*stretchFac;

  double sin_p = sin(targetLat);
  double cos_p = cos(targetLat);

  double sin_lat;
  double cos_lat;
  double lat_t;

  if ( std::abs(c2m1) > 1.0e-7 ) {
    sin_lat = sin(lonlat[LAT]);
    lat_t = asin( (c2m1+c2p1*sin_lat)/(c2p1+c2m1*sin_lat) );
  } else {         // no stretching
    lat_t = lonlat[LAT];
  }

  sin_lat = sin(lat_t);
  cos_lat = cos(lat_t);
  double sin_o = -(sin_p*sin_lat + cos_p*cos_lat*cos(lonlat[LON]));

  if ( (1.-std::abs(sin_o)) < 1.0e-7 ) {    // poles
    lonlat[LON] = 0.0;
    lonlat[LAT] = copysign( 0.5*M_PI, sin_o );
  } else {
    lonlat[LAT] = asin( sin_o );
    lonlat[LON] = targetLon + atan2( -cos_lat*sin(lonlat[LON]), -sin_lat*cos_p+cos_lat*sin_p*cos(lonlat[LON]));
    if ( lonlat[LON] < 0.0 ) {
      lonlat[LON] = lonlat[LON] + 2.0*M_PI;
    } else if ( lonlat[LON] >= 2.0*M_PI ) {
      lonlat[LON] = lonlat[LON] - 2.0*M_PI;
    }
  }

}

// -------------------------------------------------------------------------------------------------

}  // namespace detail
}  // namespace projection
}  // namespace atlas
