/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "CubedSphereProjectionBase.h"

#include <cmath>
#include <functional>
#include <array>
#include <limits>
#include <iostream>
#include <iomanip>
#include <memory>

#include "atlas/projection/detail/ProjectionUtilities.h"
#include "atlas/runtime/Log.h"
#include "atlas/runtime/Trace.h"
#include "atlas/runtime/Exception.h"
#include "atlas/util/Constants.h"
#include "atlas/util/CoordinateEnums.h"

//#include "atlas/grid/Tiles.h"

namespace atlas {
namespace projection {
namespace detail {

// -------------------------------------------------------------------------------------------------
// Helper functions and variables local to this translation unit
namespace {

static constexpr bool debug = true; // constexpr so compiler can optimize `if ( debug ) { ... }` out

static constexpr double deg2rad = util::Constants::degreesToRadians();
static constexpr double rad2deg = util::Constants::radiansToDegrees();

static bool is_tiny( const double& x ) {
    constexpr double epsilon = 1.e-15;
    return (std::abs(x) < epsilon );
}

static bool is_same( const double& x, const double& y, const double& tol = 1.0 ) {
    constexpr double epsilon = 1.e-15;
    return (std::abs(x-y) < epsilon * tol);
}

static double abs_diff( const double& x, const double& y ) {
    return (std::abs(x-y));
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

static const std::array<std::function<void(double[])>, 6> tileRotate = {
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

static const std::array<std::function<void(double[])>,6> tileRotateInverse = {
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
        sin_lat = std::sin(lonlat[LAT]);
        lat_t = std::asin( (c2m1+c2p1*sin_lat)/(c2p1+c2m1*sin_lat) );
    } else {         // no stretching
        lat_t = lonlat[LAT];
    }

    sin_lat = std::sin(lat_t);
    cos_lat = std::cos(lat_t);
    double sin_o = -(sin_p*sin_lat + cos_p*cos_lat*cos(lonlat[LON]));

    if ( (1.-std::abs(sin_o)) < 1.0e-7 ) {    // poles
        lonlat[LON] = 0.0;
        lonlat[LAT] = std::copysign( 0.5*M_PI, sin_o );
    } else {
        lonlat[LAT] = std::asin( sin_o );
        lonlat[LON] = targetLon + atan2( -cos_lat*std::sin(lonlat[LON]), -sin_lat*cos_p+cos_lat*sin_p*std::cos(lonlat[LON]));
        if ( lonlat[LON] < 0.0 ) {
            lonlat[LON] = lonlat[LON] + 2.0*M_PI;
        } else if ( lonlat[LON] >= 2.0*M_PI ) {
            lonlat[LON] = lonlat[LON] - 2.0*M_PI;
        }
    }
}

void sphericalToCartesian(const double lonlat[], double xyz[] ) {
    auto crd_sys = ProjectionUtilities::CoordinateSystem::LEFT_HAND;
    constexpr double radius = 1.;
    ProjectionUtilities::sphericalToCartesian(lonlat, xyz, crd_sys, radius);
}

void cartesianToSpherical(const double xyz[], double lonlat[] ) {
    auto crd_sys = ProjectionUtilities::CoordinateSystem::LEFT_HAND;
    constexpr double radius = 0.; // --> equivalent to radius = norm(xyz)
    ProjectionUtilities::cartesianToSpherical(xyz, lonlat, crd_sys, radius);
}

} // namespace


// -------------------------------------------------------------------------------------------------

CubedSphereProjectionBase::CubedSphereProjectionBase( const eckit::Parametrisation& params ) :
    CubedSphereTiles_(params)
{
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

void CubedSphereProjectionBase::xy2lonlat_post( double xyz[], const idx_t & t, double crd[] ) const {

    cartesianToSpherical(xyz, crd);

    if (crd[LON] < 0.0) { crd[LON] += 2.0*M_PI; }
    crd[LON] = crd[LON] - M_PI;

    if ( debug ) {
        Log::debug() << "xy2lonlat:: lonlat before rotation : "  << crd[LON] << " " << crd[LAT]  << std::endl;
    }

    // Convert to cartesian
    sphericalToCartesian(crd, xyz);

    // Perform tile specific rotation
    tileRotate[t](xyz);

    // Back to lonlat
    cartesianToSpherical(xyz, crd);

    // Shift longitude
    if (shiftLon_ != 0.0) {
      crd[LON] = crd[LON] + shiftLon_*deg2rad;
      if (crd[LON] < -M_PI) {crd[LON] =  2*M_PI + crd[LON];}
      if (crd[LON] >  M_PI) {crd[LON] = -2*M_PI + crd[LON];}
    }

    // To 0, 360
    if (crd[LON] < 0.0) { crd[LON] = 2.*M_PI + crd[LON]; }

    // Schmidt transform
    if (doSchmidt_) {
       schmidtTransform(stretchFac_,
                        targetLon_*deg2rad,
                        targetLat_*deg2rad,
                        crd);
    }

    // longitude does not make sense at the poles - set to 0.
    if (std::abs(std::abs(crd[LAT]) - M_PI_2) < 1e-15) { crd[LON] = 0.; }

    crd[LON] *= rad2deg;
    crd[LAT] *= rad2deg;
}

// -------------------------------------------------------------------------------------------------

void CubedSphereProjectionBase::lonlat2xy_pre( double crd[], idx_t & t, double xyz[] ) const {

    if (std::abs(crd[LON]) < 1e-15) { crd[LON] = 0.; }
    if (std::abs(crd[LAT]) < 1e-15) { crd[LAT] = 0.; }

    // convert degrees to radians
    crd[LON] *= deg2rad;
    crd[LAT] *= deg2rad;

    // To [-pi/4, 7/4 * pi)
    if (crd[LON] >= 1.75 * M_PI) { crd[LON] += -2.*M_PI; }

    // find tile which this lonlat is linked to
    // works [-pi/4, 7/4 * pi)
    t = CubedSphereTiles_.tileFromLonLat(crd);

    sphericalToCartesian(crd, xyz);
    tileRotateInverse[t](xyz);

}

// -------------------------------------------------------------------------------------------------


void CubedSphereProjectionBase::xy2alphabetat(const double xy[], idx_t& t, double ab[]) const {
    // xy is in degrees while ab is in radians
    // ab are the  (alpha, beta) coordinates and t is the tile index.
    static std::array<double,6> xOffset{0., 1., 1., 2., 3., 3.}; // could become constexpr with C++17
    static std::array<double,6> yOffset{1., 1., 2., 1., 1., 0.};

    t = CubedSphereTiles_.tileFromXY(xy);
    double normalisedX = xy[XX]/90.;
    double normalisedY = (xy[YY] + 135.)/90.;
    ab[LON] = (normalisedX - xOffset[t])* M_PI_2 - M_PI_4;
    ab[LAT] = (normalisedY - yOffset[t])* M_PI_2 - M_PI_4;
}

// -------------------------------------------------------------------------------------------------

void CubedSphereProjectionBase::alphabetat2xy(const idx_t& t, const double ab[], double xy[]) const {
    // xy is in degrees while ab is in radians
    // (alpha, beta) and tiles.
    static std::array<double,6> xOffset{0., 90., 90., 180., 270., 270.}; // could become constexpr with C++17
    static std::array<double,6> yOffset{-45., -45., 45., -45., -45., -135.};
    double normalisedX = (ab[LON] + M_PI_4)/M_PI_2;
    double normalisedY = (ab[LAT] + M_PI_4)/M_PI_2;
    xy[XX] = normalisedX * 90. + xOffset[t];
    xy[YY] = normalisedY * 90. + yOffset[t];

    CubedSphereTiles_.enforceXYdomain(xy);
}



// -------------------------------------------------------------------------------------------------

}  // namespace detail
}  // namespace projection
}  // namespace atlas
