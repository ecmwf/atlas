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

  // -----------------------------------------------------------------------------------------------

  static void cartesianToSpherical(double xyz[], double lonlat[], bool right_hand = false) {
    // Left or right hand system
    double zz_fac = 1.0;
    if (right_hand) {zz_fac = 0.0;}

    double r = sqrt(xyz[XX]*xyz[XX] + xyz[YY]*xyz[YY] + xyz[ZZ]*xyz[ZZ]);
    lonlat[LON] = atan2(xyz[YY],xyz[XX]);
    if ( (abs(xyz[XX]) + abs(xyz[YY])) < util::Constants::esl() ) {
      lonlat[LON] = 0.0;  // Poles:
    }
    lonlat[LAT] = acos(xyz[ZZ]/r) - zz_fac * M_PI/2.0;
  }

  //------------------------------------------------------------------------------------------------

  static void sphericalToCartesian(double lonlat[], double xyz[], bool right_hand = false,
                                  bool unit = false) {
    // Radius
    double r = util::Earth::radius();
    if (unit) {r = 1.0;}

    // Left or right hand system
    double zz_fac = 1.0;
    if (right_hand) {zz_fac = -1.0;}

    xyz[XX] = r * cos(lonlat[LON]) * cos(lonlat[LAT]);
    xyz[YY] = r * sin(lonlat[LON]) * cos(lonlat[LAT]);
    xyz[ZZ] = zz_fac * -r * sin(lonlat[LAT]);
  }

  //------------------------------------------------------------------------------------------------

  static void rotate3dX(double angle, double xyz[]) {
    const double c = cos(angle);
    const double s = sin(angle);
    double xyz_in[3];
    std::copy(xyz, xyz+3, xyz_in);
    xyz[YY] =  c*xyz_in[YY] + s*xyz_in[ZZ];
    xyz[ZZ] = -s*xyz_in[YY] + c*xyz_in[ZZ];
  };

  //------------------------------------------------------------------------------------------------

  static void rotate3dY(double angle, double xyz[]) {
    const double c = cos(angle);
    const double s = sin(angle);
    double xyz_in[3];
    std::copy(xyz, xyz+3, xyz_in);
    xyz[XX] = c*xyz_in[XX] - s*xyz_in[ZZ];
    xyz[ZZ] = s*xyz_in[XX] + c*xyz_in[ZZ];
  };

  //------------------------------------------------------------------------------------------------

  static void rotate3dZ(double angle, double xyz[]) {
    const double c = cos(angle);
    const double s = sin(angle);
    double xyz_in[3];
    std::copy(xyz, xyz+3, xyz_in);
    xyz[XX] =  c*xyz_in[XX] + s*xyz_in[YY];
    xyz[YY] = -s*xyz_in[XX] + c*xyz_in[YY];
  };

  //------------------------------------------------------------------------------------------------

  static void vect_cross(double xyz1[], double xyz2[], double xyzp[]) {
    // Perform cross products of 3D vectors: xyzp = xyz1 X xyz2
    xyzp[XX] = xyz1[YY] * xyz2[ZZ] - xyz1[ZZ] * xyz2[YY];
    xyzp[YY] = xyz1[ZZ] * xyz2[XX] - xyz1[XX] * xyz2[ZZ];
    xyzp[ZZ] = xyz1[XX] * xyz2[YY] - xyz1[YY] * xyz2[XX];
  }

  //------------------------------------------------------------------------------------------------

  static void mirror_latlon(double lonlat2[], double lonlat3[], double lonlat1[], double lonlat4[]) {
    // Given the "mirror" as defined by (lon1, lat1), (lon2, lat2), and center
    // of the sphere, compute the mirror image of (lon0, lat0) as  (lon3, lat3)

    double xyz0[3];
    double xyz1[3];
    double xyz2[3];

    sphericalToCartesian(lonlat1, xyz0, true, true);
    sphericalToCartesian(lonlat2, xyz1, true, true);
    sphericalToCartesian(lonlat3, xyz2, true, true);

    double xyz1xyz2[3];
    vect_cross(xyz1, xyz2, xyz1xyz2);

    double dot = sqrt(xyz1xyz2[XX]*xyz1xyz2[XX] +
                      xyz1xyz2[YY]*xyz1xyz2[YY] +
                      xyz1xyz2[ZZ]*xyz1xyz2[ZZ]);

    for (int n = 0; n < 3; n++) {
      xyz1xyz2[n] = xyz1xyz2[n] / dot;
    }

    dot = xyz1xyz2[XX]*xyz0[XX] + xyz1xyz2[YY]*xyz0[YY] + xyz1xyz2[ZZ]*xyz0[ZZ];

    double xyzout[3];

    for (int n = 0; n < 3; n++) {
      xyzout[n] = xyz0[n] - 2.0*dot*xyz1xyz2[n];
    }

    cartesianToLatLon(xyzout[0], xyzout[1], xyzout[2], lonlat4);
  }

  //--------------------------------------------------------------------------------------------------

};

//--------------------------------------------------------------------------------------------------

} // detail
} // projection
} // atlas
