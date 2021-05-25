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
#include "atlas/runtime/Log.h"
#include "atlas/util/Constants.h"
#include "atlas/util/CoordinateEnums.h"

namespace atlas {
namespace projection {
namespace detail {

// -------------------------------------------------------------------------------------------------
// Helper functions and variables local to this translation unit
namespace {
static constexpr bool debug = false; // constexpr so compiler can optimize `if ( debug ) { ... }` out

static bool is_tiny( const double& x ) {
    constexpr double epsilon = 1.e-15;
    return (std::abs(x) < epsilon );
}

static bool is_same( const double& x, const double& y ) {
    constexpr double epsilon = 1.e-15;
    return (std::abs(x-y) < epsilon );
}

// --- Functions for xy to latlon on each tile

static void tile1Rotate( double xyz[] ) {
    //  Face 1, no rotation.
}
static void tile2Rotate( double xyz[] ) {
    //  Face 2: rotate -90.0 degrees about z axis
    constexpr double angle = -M_PI / 2.0;
    ProjectionUtilities::rotate3dZ(angle, xyz);
}
static void tile3Rotate( double xyz[] ) {
    //  Face 3: rotate -90.0 degrees about z axis
    //          rotate  90.0 degrees about x axis
    constexpr double angle_z = -M_PI / 2.0;
    ProjectionUtilities::rotate3dZ(angle_z, xyz);
    constexpr double angle_x = M_PI / 2.0;
    ProjectionUtilities::rotate3dX(angle_x, xyz);
}
static void tile4Rotate( double xyz[] ) {
    //  Face 4: rotate -180.0 degrees about z axis
    constexpr double angle = -M_PI;
    ProjectionUtilities::rotate3dZ(angle, xyz);
}
static void tile5Rotate( double xyz[] ) {
    //  Face 5: rotate 90.0 degrees about z axis
    constexpr double angle = M_PI / 2.0;
    ProjectionUtilities::rotate3dZ(angle, xyz);
}
static void tile6Rotate( double xyz[] ) {
    // Face 6: rotate 90.0 degrees about y axis
    //         rotate 90.0 degrees about z axis
    constexpr double angle_y = M_PI / 2.0;
    ProjectionUtilities::rotate3dY(angle_y, xyz);
    constexpr double angle_z = M_PI / 2.0;
    ProjectionUtilities::rotate3dZ(angle_z, xyz);
}

static const std::vector<std::function<void(double[])>> tileRotate = {
    [](double xyz[]){tile1Rotate(xyz);},
    [](double xyz[]){tile2Rotate(xyz);},
    [](double xyz[]){tile3Rotate(xyz);},
    [](double xyz[]){tile4Rotate(xyz);},
    [](double xyz[]){tile5Rotate(xyz);},
    [](double xyz[]){tile6Rotate(xyz);}
};

// Functions for latlon to xy on each tile
static void tile1RotateInverse( double xyz[] ) {
    //  Face 1, no rotation.
}
static void tile2RotateInverse( double xyz[] ) {
    //  Face 2: rotate 90.0 degrees about z axis
    constexpr double angle = M_PI / 2.0;
    ProjectionUtilities::rotate3dZ(angle, xyz);
}
static void tile3RotateInverse( double xyz[] ) {
    //  Face 3: rotate  90.0 degrees about x axis
    //          rotate -90.0 degrees about z axis
    constexpr double angle_x = -M_PI / 2.0;
    ProjectionUtilities::rotate3dX(angle_x, xyz);
    constexpr double angle_z = M_PI / 2.0;
    ProjectionUtilities::rotate3dZ(angle_z, xyz);
}
static void tile4RotateInverse( double xyz[] ) {
    //  Face 4: rotate  180.0 degrees about z axis
    constexpr double angle = M_PI;
    ProjectionUtilities::rotate3dZ(angle, xyz);
}
static void tile5RotateInverse( double xyz[] ) {
    //  Face 5: rotate -90.0 degrees about y axis
    constexpr double angle = - M_PI / 2.0;
    ProjectionUtilities::rotate3dZ(angle, xyz);
}
static void tile6RotateInverse( double xyz[] ) {
    //  Face 6: rotate -90.0 degrees about y axis
    //          rotate -90.0 degrees about z axis
    constexpr double angle_z = -M_PI / 2.;
    ProjectionUtilities::rotate3dZ(angle_z, xyz);
    constexpr double angle_y = -M_PI / 2.;
    ProjectionUtilities::rotate3dY(angle_y, xyz);
}

static const std::vector<std::function<void(double[])>> tileRotateInverse = {
    [](double xyz[]){tile1RotateInverse(xyz);},
    [](double xyz[]){tile2RotateInverse(xyz);},
    [](double xyz[]){tile3RotateInverse(xyz);},
    [](double xyz[]){tile4RotateInverse(xyz);},
    [](double xyz[]){tile5RotateInverse(xyz);},
    [](double xyz[]){tile6RotateInverse(xyz);},
};

static void schmidtTransform( double stretchFac, double targetLon,
                              double targetLat, double lonlat[]) {


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


} // namespace


// -------------------------------------------------------------------------------------------------

CubedSphereProjectionBase::CubedSphereProjectionBase( const eckit::Parametrisation& params ) {
  ATLAS_TRACE( "CubedSphereProjectionBase::CubedSphereProjectionBase" );

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

    if (crd[LON] < 0.0) { crd[LON] += 2.0*M_PI; }
    crd[LON] = crd[LON] - M_PI;

    if ( debug ) {
        Log::debug() << "xy2lonlat:: lonlat before rotation : "  << crd[LON] << " " << crd[LAT]  << std::endl;
    }

    // Convert to cartesian
    ProjectionUtilities::sphericalToCartesian(crd, xyz, false, true);

    // Perform tile specific rotation
    tileRotate[t](xyz);

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
       schmidtTransform(stretchFac_,
                        targetLon_*atlas::util::Constants::degreesToRadians(),
                        targetLat_*atlas::util::Constants::degreesToRadians(),
                        crd);
    }

    // longitude does not make sense at the poles - set to 0.
    if (std::abs(std::abs(crd[LAT]) - M_PI_2) < 1e-15) { crd[LON] = 0.; }

    crd[LON] *= Constants::radiansToDegrees();
    crd[LAT] *= Constants::radiansToDegrees();
}

// -------------------------------------------------------------------------------------------------

void CubedSphereProjectionBase::lonlat2xypre( double crd[], idx_t & t, double xyz[] ) const {

    using util::Constants;

    if (std::abs(crd[LON]) < 1e-15) { crd[LON] = 0.; }
    if (std::abs(crd[LAT]) < 1e-15) { crd[LAT] = 0.; }

    // convert degrees to radians
    crd[LON] *= Constants::degreesToRadians();
    crd[LAT] *= Constants::degreesToRadians();

    // To [-pi/4, 7/4 * pi)
    if (crd[LON] >= 1.75 * M_PI) { crd[LON] += -2.*M_PI; }

    // find tile which this lonlat is linked to
    // works [-pi/4, 7/4 * pi)
    t = CubedSphereProjectionBase::tileFromLonLat(crd);

    ProjectionUtilities::sphericalToCartesian(crd, xyz, false, true);
    tileRotateInverse[t](xyz);

}

// -------------------------------------------------------------------------------------------------

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


  if ((xy[XX] >= 0.) && ( xy[YY] >= -45.) && (xy[XX] < 90.) && (xy[YY] < 45.)) {
     t = 0;
  } else if ((xy[XX] >= 90.) && ( xy[YY] >= -45.) && (xy[XX] < 180.) && (xy[YY] < 45.)) {
     t = 1;
  } else if ((xy[XX] >= 90.) && ( xy[YY] >= 45.) && (xy[XX] < 180.) && (xy[YY] < 135.)) {
     t = 2;
  } else if ((xy[XX] >= 180.) && ( xy[YY] > -45.) && (xy[XX] < 270.) && (xy[YY] <= 45.)) {
     t = 3;
  } else if ((xy[XX] >= 270.) && ( xy[YY] > -45.) && (xy[XX] < 360.) && (xy[YY] <= 45.)) {
     t = 4;
  } else if ((xy[XX] >= 270.) && ( xy[YY] > -135.) && (xy[XX] < 360.) && (xy[YY] <= -45.)) {
     t = 5;
  }

  // extra points
  if ( is_same( xy[XX], 0.   ) && is_same( xy[YY],  45. ) ) { t = 0; }
  if ( is_same( xy[XX], 180. ) && is_same( xy[YY], -45. ) ) { t = 1; }

  // for end iterator !!!!
  if ( is_same( xy[XX], 360. ) && is_same( xy[YY], -135.) ) { t = 5; }

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

    double zPlusAbsX = xyz[ZZ] + std::abs(xyz[XX]);
    double zPlusAbsY = xyz[ZZ] + std::abs(xyz[YY]);
    double zMinusAbsX = xyz[ZZ] - std::abs(xyz[XX]);
    double zMinusAbsY = xyz[ZZ] - std::abs(xyz[YY]);

    // Note that this method can lead to roundoff errors that can
    // cause the tile selection to fail.
    // To this end we enforce that tiny values close (in roundoff terms)
    // to a boundary should end up exactly on the boundary.
    if ( is_tiny(zPlusAbsX) ) { zPlusAbsX = 0.; }
    if ( is_tiny(zPlusAbsY) ) { zPlusAbsY = 0.; }
    if ( is_tiny(zMinusAbsX) ) { zMinusAbsX = 0.; }
    if ( is_tiny(zMinusAbsY) ) { zMinusAbsY = 0.; }

    if (lon >= 1.75 * M_PI  || lon < 0.25 * M_PI) {
        if  ( (zPlusAbsX <= 0.) && (zPlusAbsY <= 0.) ) {
           t = 2;
        } else if ( (zMinusAbsX > 0.) && (zMinusAbsY > 0.) ) {
           t = 5;
        } else {
           t = 0;
        }
        // extra point corner point
        if ( is_same( lon, -0.25 * M_PI ) && is_same( lat, cornerLat ) ) { t = 0; }
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
        if ( is_same(lon, 0.75 * M_PI) && is_same( lat, -cornerLat ) ) { t = 1; }
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

    if( debug ) {
        Log::debug() << "tileFromLonLat:: lonlat abs xyz t = "
                     << std::setprecision(std::numeric_limits<double>::digits10 + 1)
                     << zPlusAbsX << " " << zPlusAbsY << " " << zMinusAbsX << " " << zMinusAbsY << " "
                     << crd[LON] << " " << crd[LAT] << " "
                     << xyz[XX] << " " << xyz[YY] << " " << xyz[ZZ] << " " << t
                     << std::endl;
    }

    return t;
}

// -------------------------------------------------------------------------------------------------

}  // namespace detail
}  // namespace projection
}  // namespace atlas
