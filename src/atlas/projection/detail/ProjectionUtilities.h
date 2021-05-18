/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <cmath>
#include <functional>

#include "atlas/domain.h"
#include "atlas/util/Constants.h"
#include "atlas/util/CoordinateEnums.h"
#include "atlas/util/Earth.h"
#include "atlas/util/Point.h"
#include "eckit/geometry/Sphere.h"

namespace atlas {
namespace projection {
namespace detail {

// -------------------------------------------------------------------------------------------------

struct ProjectionUtilities {

  // -----------------------------------------------------------------------------------------------

  static constexpr double epsilon() { return 1.0e-10; }

  // -----------------------------------------------------------------------------------------------

  static void cartesianToSpherical(const double xyz[], double lonlat[],
                                   const bool right_hand = false) {

    using eckit::geometry::Sphere;
    using util::Constants;

    // Make point objects.
    const auto pointXyz = PointXYZ(xyz);
    auto pointLonLat = PointLonLat();

    // Transform coordinates.
    Sphere::convertCartesianToSpherical(
      PointXYZ::norm(pointXyz), pointXyz, pointLonLat);

    // Copy to array.
    lonlat[LON] = pointLonLat.lon() * Constants::degreesToRadians();
    lonlat[LAT] = -pointLonLat.lat() * Constants::degreesToRadians();

    // Left or right hand system.
    if (right_hand) lonlat[LAT] += 0.5 * M_PI;

  }

  //------------------------------------------------------------------------------------------------

  static void sphericalToCartesian(const double lonlat[], double xyz[],
                                   const bool right_hand = false, const bool unit = false) {

    using eckit::geometry::Sphere;
    using util::Constants;

    // Make point objects.
    const auto pointLonLat = PointLonLat(lonlat) * Constants::radiansToDegrees();
    auto pointXyz = PointXYZ();

    // Set Radius
    auto r = util::Earth::radius();
    if (unit) {r = 1.0;}

    // Transform coordinates.
    Sphere::convertSphericalToCartesian(r, pointLonLat, pointXyz, 0.0);

    // Copy to array.
    xyz[XX] = pointXyz.x();
    xyz[YY] = pointXyz.y();
    xyz[ZZ] = pointXyz.z();

    // Left or right hand system
    if (!right_hand) xyz[ZZ] *= -1;

  }

  //------------------------------------------------------------------------------------------------

  static void rotate3dX(const double angle, double xyz[]) {
    const double c = std::cos(angle);
    const double s = std::sin(angle);
    double xyz_in[3];
    std::copy(xyz, xyz+3, xyz_in);
    xyz[YY] =  c*xyz_in[YY] + s*xyz_in[ZZ];
    xyz[ZZ] = -s*xyz_in[YY] + c*xyz_in[ZZ];
  };

  //------------------------------------------------------------------------------------------------

  static void rotate3dY(const double angle, double xyz[]) {
    const double c = std::cos(angle);
    const double s = std::sin(angle);
    double xyz_in[3];
    std::copy(xyz, xyz+3, xyz_in);
    xyz[XX] = c*xyz_in[XX] - s*xyz_in[ZZ];
    xyz[ZZ] = s*xyz_in[XX] + c*xyz_in[ZZ];
  };

  //------------------------------------------------------------------------------------------------

  static void rotate3dZ(const double angle, double xyz[]) {
    const double c = std::cos(angle);
    const double s = std::sin(angle);
    double xyz_in[3];
    std::copy(xyz, xyz+3, xyz_in);
    xyz[XX] =  c*xyz_in[XX] + s*xyz_in[YY];
    xyz[YY] = -s*xyz_in[XX] + c*xyz_in[YY];
  };

  //------------------------------------------------------------------------------------------------

};

//--------------------------------------------------------------------------------------------------

} // detail
} // projection
} // atlas
