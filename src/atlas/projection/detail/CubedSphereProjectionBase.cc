/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <cmath>
#include <limits>
#include <iostream>
#include <iomanip>
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

void CubedSphereProjectionBase::xy2lonlatpost( double xyz[], const idx_t & t, double crd[] ) const {

    using util::Constants;

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
    if (shiftLon_ != 0.0) {
      crd[LON] = crd[LON] + shiftLon_*atlas::util::Constants::degreesToRadians();
      if (crd[LON] < -M_PI) {crd[LON] =  2*M_PI + crd[LON];}
      if (crd[LON] >  M_PI) {crd[LON] = -2*M_PI + crd[LON];}
    }

    // To 0, 360
    if (crd[LON] < 0.0) { crd[LON] = 2.*M_PI + crd[LON]; }

    // Schmidt transform
    if (doSchmidt_) {
       this->schmidtTransform(stretchFac_,
                              targetLon_*atlas::util::Constants::degreesToRadians(),
                              targetLat_*atlas::util::Constants::degreesToRadians(),
                              crd);
    }

    // longitude does not make sense at the poles - set to 0.
    if (std::abs(std::abs(crd[LAT]) - M_PI_2) < 1e-15) crd[LON] = 0.;

    crd[LON] *= Constants::radiansToDegrees();
    crd[LAT] *= Constants::radiansToDegrees();
}

// -------------------------------------------------------------------------------------------------

void CubedSphereProjectionBase::lonlat2xypre( double crd[], idx_t & t, double xyz[] ) const {

    using util::Constants;

    if (std::abs(crd[LON]) < 1e-15) crd[LON] = 0.;
    if (std::abs(crd[LAT]) < 1e-15) crd[LAT] = 0.;

    // convert degrees to radians
    crd[LON] *= Constants::degreesToRadians();
    crd[LAT] *= Constants::degreesToRadians();

    // To [-pi/4, 7/4 * pi)
    if (crd[LON] >= 1.75 * M_PI) crd[LON] += -2.*M_PI;

    // find tile which this lonlat is linked to
    // works [-pi/4, 7/4 * pi)
    t = CubedSphereProjectionBase::tileFromLonLat(crd);

    ProjectionUtilities::sphericalToCartesian(crd, xyz, false, true);
    tileRotateInverse.at(t)(xyz);

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
  double angle;
  angle = -M_PI;
  ProjectionUtilities::rotate3dZ(angle, xyz);
  std::cout << "tile4 Rotate after 3dZ" << xyz[0]  << " " << xyz[1] << " " << xyz[2] << std::endl;
}

// -------------------------------------------------------------------------------------------------

void CubedSphereProjectionBase::tile5Rotate( double xyz[] ) const {
  //  Face 5: rotate 90.0 degrees about z axis
  double angle;
  angle = M_PI / 2.0;
  ProjectionUtilities::rotate3dZ(angle, xyz);
}

// -------------------------------------------------------------------------------------------------

void CubedSphereProjectionBase::tile6Rotate( double xyz[] ) const {
  // Face 6: rotate 90.0 degrees about y axis
  //         rotate 90.0 degrees about z axis
  double angle;
  angle = M_PI / 2.0;
  ProjectionUtilities::rotate3dY(angle, xyz);
  angle = M_PI / 2.0;
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
  //  Face 4: rotate  180.0 degrees about z axis
  double angle;
  angle = M_PI;
  ProjectionUtilities::rotate3dZ(angle, xyz);
}

// -------------------------------------------------------------------------------------------------

void CubedSphereProjectionBase::tile5RotateInverse( double xyz[] ) const {
  //  Face 5: rotate -90.0 degrees about y axis
  double angle;
  angle = - M_PI / 2.0;
  ProjectionUtilities::rotate3dZ(angle, xyz);
}

// -------------------------------------------------------------------------------------------------

void CubedSphereProjectionBase::tile6RotateInverse( double xyz[] ) const {
  //  Face 6: rotate -90.0 degrees about y axis
  //          rotate -90.0 degrees about z axis
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

  double sin_p = std::sin(targetLat);
  double cos_p = std::cos(targetLat);

  double sin_lat;
  double cos_lat;
  double lat_t;

  if ( std::abs(c2m1) > 1.0e-7 ) {
    sin_lat = sin(lonlat[LAT]);
    lat_t = std::asin( (c2m1+c2p1*sin_lat)/(c2p1+c2m1*sin_lat) );
  } else {         // no stretching
    lat_t = lonlat[LAT];
  }

  sin_lat = std::sin(lat_t);
  cos_lat = std::cos(lat_t);
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

idx_t CubedSphereProjectionBase::tileFromXY( const double xy[] ) const {
  // Assume one face-edge is of length 90 degrees.
  //
  //   y ^
  //     |
  //    135              ----------
  //     |              |     ^    |
  //     |              |          |
  //     |              |=<   2   <|
  //     |              |     v    |
  //     |              |     =    |
  //     45  0----------2----------3----------4----------
  //     |   |    ^     |     ^    |    =     |     =    |
  //     |   |          |          |    ^     |     ^    |
  //     |   |=<  0    <|=<   1   <|=<  3    <|=<   4   <|
  //     |   |    v     |     v    |          |          |
  //     |   |    =     |     =    |    v     |     v    |
  //    -45  0 ---------1----------1----------5----------
  //     |                                    |     =    |
  //     |                                    |     ^    |
  //     |                                    |=<   5   <|
  //     |                                    |          |
  //     |                                    |     v    |
  //   -135                                    ----------(5 for end iterator)
  //     ----0---------90--------180--------270--------360--->  x

  idx_t t{-1};

  if ((xy[0] >= 0.) && ( xy[1] >= -45.) && (xy[0] < 90.) && (xy[1] < 45.)) {
     t = 0;
  } else if ((xy[0] >= 90.) && ( xy[1] >= -45.) && (xy[0] < 180.) && (xy[1] < 45.)) {
     t = 1;
  } else if ((xy[0] >= 90.) && ( xy[1] >= 45.) && (xy[0] < 180.) && (xy[1] < 135.)) {
     t = 2;
  } else if ((xy[0] >= 180.) && ( xy[1] > -45.) && (xy[0] < 270.) && (xy[1] <= 45.)) {
     t = 3;
  } else if ((xy[0] >= 270.) && ( xy[1] > -45.) && (xy[0] < 360.) && (xy[1] <= 45.)) {
     t = 4;
  } else if ((xy[0] >= 270.) && ( xy[1] > -135.) && (xy[0] < 360.) && (xy[1] <= -45.)) {
     t = 5;
  }

  // extra points
  if ((std::abs(xy[0]) < 1e-13) && (std::abs(xy[1] - 45.) < 1e-13)) t = 0;
  if ((std::abs(xy[0] - 180.) < 1e-13) && (std::abs(xy[1] + 45.) < 1e-13)) t = 1;

  // for end iterator !!!!
  if ((std::abs(xy[0] - 360.) < 1e-13) && (std::abs(xy[1] + 135.) < 1e-13)) t = 5;

  return t;
}

// assuming that crd[] here holds the latitude/longitude in RADIANS.
// expecting longitude range between 0 and 2* PI
idx_t CubedSphereProjectionBase::tileFromLonLat(const double crd[]) const {
    idx_t t(-1);
    double xyz[3];

    ProjectionUtilities::sphericalToCartesian(crd, xyz, false, true);

    const double cornerLat= std::asin(1./std::sqrt(3.0));
      // the magnitude of the latitude at the corners of the cube (not including the sign)
      // in radians.

    const double & lon = crd[LON];
    const double & lat = crd[LAT];

    double zPlusAbsX = xyz[2] + abs(xyz[0]);
    double zPlusAbsY = xyz[2] + abs(xyz[1]);
    double zMinusAbsX = xyz[2] - abs(xyz[0]);
    double zMinusAbsY = xyz[2] - abs(xyz[1]);

    // Note that this method can lead to roundoff errors that can
    // cause the tile selection to fail.
    // To this end we enforce that tiny values close (in roundoff terms)
    // to a boundary should end up exactly on the boundary.
    double tolerance = 1e-15;
    if (abs(zPlusAbsX) < tolerance) zPlusAbsX = 0.;
    if (abs(zPlusAbsY) < tolerance) zPlusAbsY = 0.;
    if (abs(zMinusAbsX) < tolerance) zMinusAbsX = 0.;
    if (abs(zMinusAbsY) < tolerance) zMinusAbsY = 0.;

    if (lon >= 1.75 * M_PI  || lon < 0.25 * M_PI) {
        if  ( (zPlusAbsX <= 0.) && (zPlusAbsY <= 0.) ) {
           t = 2;
        } else if ( (zMinusAbsX > 0.) && (zMinusAbsY > 0.) ) {
           t = 5;
        } else {
           t = 0;
        }
        // extra point corner point
        if ( (abs(lon + 0.25 * M_PI) < 1e-13) &&
             (abs(lat - cornerLat) < 1e-13) ) t = 0;
    }

    if (lon >= 0.25 * M_PI  && lon < 0.75 * M_PI) {
        // interior
        if  ( (zPlusAbsX <= 0.) && (zPlusAbsY <= 0.) ) {
            t = 2;
        } else if  ( (zMinusAbsX > 0.) && (zMinusAbsY > 0.) ) {
            t = 5;
        } else {
            t = 1;
        }
    }

    if (lon >= 0.75 * M_PI  && lon < 1.25 * M_PI) {
        // interior
        if  ( (zPlusAbsX < 0.) && (zPlusAbsY < 0.) ) {
            t = 2;
        } else if  ( (zMinusAbsX >= 0.) && (zMinusAbsY >= 0.) ) {
            t = 5;
        } else {
            t = 3;
        }
        // extra point corner point
        if ( (abs(lon - 0.75 * M_PI) < 1e-13) &&
             (abs(lat + cornerLat) < 1e-13) ) t = 1;
    }

    if (lon >= 1.25 * M_PI  && lon < 1.75 * M_PI) {
        // interior
        if  ( (zPlusAbsX < 0.) && (zPlusAbsY < 0.) ) {
            t = 2;
        } else if  ( (zMinusAbsX >= 0.) && (zMinusAbsY >= 0.) ) {
            t = 5;
        } else {
            t = 4;
        }
    }

    std::cout << "tileFromLonLat:: lonlat abs xyz t = " <<
                 std::setprecision(std::numeric_limits<double>::digits10 + 1) <<
        zPlusAbsX  << " " << zPlusAbsY << " " << zMinusAbsX << " " << zMinusAbsY << " " <<
        crd[0] << " " << crd[1]  << " "  <<
        xyz[0] << " " << xyz[1]  << " " << xyz[2] << " " << t << std::endl;

    return t;
}

// -------------------------------------------------------------------------------------------------

}  // namespace detail
}  // namespace projection
}  // namespace atlas
