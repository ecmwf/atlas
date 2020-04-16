/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <functional>

#include "atlas/domain.h"
#include "atlas/util/Constants.h"
#include "atlas/util/CoordinateEnums.h"
#include "atlas/util/Earth.h"

namespace atlas {
namespace projection {
namespace detail {

// -------------------------------------------------------------------------------------------------

struct ProjectionUtilities {

  static void cartesianToLatLon(double p1, double p2, double p3, double lonlat[]) {
    // Input is cartesian coordinates, output is degrees in radians

    double dist = sqrt(pow(p1, 2) + pow(p2, 2) + pow(p3, 2));

    double lon;
    double lat;

    p1 /= dist;
    p2 /= dist;
    p3 /= dist;

    if ( (abs(p1) + abs(p2)) < util::Constants::esl() ) {
      lon = 0.0;
    } else {
      lon = atan2( p2, p1 );
    }

    // 0 to 2pi
    if (lon < 0.0) {
      lon += 2.0*M_PI;
    }

    lat = asin(p3);

    lonlat[LON] = lon;
    lonlat[LAT] = lat;
  }

  // -------------------------------------------------------------------------------------------------

  static void cartesianToSpherical(double xyz[], double lonlat[]) {
    double r = sqrt(xyz[XX]*xyz[XX] + xyz[YY]*xyz[YY] + xyz[ZZ]*xyz[ZZ]);
    lonlat[LON] = atan2(xyz[YY],xyz[XX]);
    if ( (abs(xyz[XX]) + abs(xyz[YY])) < util::Constants::esl() ) {
      lonlat[LON] = 0.0;  // Poles:
    }
    lonlat[LAT] = acos(xyz[ZZ]/r) - M_PI/2.0;
  }

  //--------------------------------------------------------------------------------------------------

  static void sphericalToCartesian(double lonlat[], double xyz[]) {
    double r = util::Earth::radius();
    xyz[XX] = r * cos(lonlat[LON]) * cos(lonlat[LAT]);
    xyz[YY] = r * sin(lonlat[LON]) * cos(lonlat[LAT]);
    xyz[ZZ] = -r * sin(lonlat[LAT]);
  }

  //--------------------------------------------------------------------------------------------------

  static void rotate3dX(double angle, double xyz[]) {
    double c = cos(angle);
    double s = sin(angle);
    auto xyz_in = xyz;
    xyz[YY] =  c*xyz_in[YY] + s*xyz_in[ZZ];
    xyz[ZZ] = -s*xyz_in[YY] + c*xyz_in[ZZ];
  };

  //--------------------------------------------------------------------------------------------------

  static void rotate3dY(double angle, double xyz[]) {
    double c = cos(angle);
    double s = sin(angle);
    auto xyz_in = xyz;
    xyz[XX] = c*xyz_in[XX] - s*xyz_in[ZZ];
    xyz[ZZ] = s*xyz_in[XX] + c*xyz_in[ZZ];
  };

  //--------------------------------------------------------------------------------------------------

  static void rotate3dZ(double angle, double xyz[]) {
    double c = cos(angle);
    double s = sin(angle);
    auto xyz_in = xyz;
    xyz[XX] =  c*xyz_in[XX] + s*xyz_in[YY];
    xyz[YY] = -s*xyz_in[XX] + c*xyz_in[YY];
  };

  //--------------------------------------------------------------------------------------------------

  // Rotate 3D vector in catesian space, angle in radians
  static std::vector<std::function<void(double, double[])>> rotate3d =
    {[](double angle, double xyz[]){rotate3dX(angle, xyz);},
     [](double angle, double xyz[]){rotate3dY(angle, xyz);},
     [](double angle, double xyz[]){rotate3dZ(angle, xyz);}
    };

};

//--------------------------------------------------------------------------------------------------

} // detail
} // projection
} // atlas
