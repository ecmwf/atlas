/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#pragma once

#include "atlas/util/Constants.h"

namespace atlas {
class PointLonLat;
class PointXYZ;
}  // namespace atlas

//------------------------------------------------------------------------------------------------------

namespace atlas {
namespace util {

//------------------------------------------------------------------------------------------------------

struct Sphere {
    // Great-circle central angle between two points, in radians
    static double centralAngle( const PointLonLat&, const PointLonLat& );
    static double centralAngle( const PointXYZ&, const PointXYZ&, const double& radius );

    // Great-circle distance between two points
    static double distanceInMeters( const PointLonLat&, const PointLonLat&, const double& radius );
    static double distanceInMeters( const PointXYZ&, const PointXYZ&, const double& radius );

    // Great-circle intermediate position provided two circle points and
    // longitude, in degrees
    static void greatCircleLatitudeGivenLongitude( const PointLonLat&, const PointLonLat&, PointLonLat& );

    // Convert spherical coordinates to Cartesian
    static PointXYZ convertSphericalToCartesian( const PointLonLat&, const double& radius, const double& height );

    // Convert Cartesian coordinates to spherical
    static PointLonLat convertCartesianToSpherical( const PointXYZ&, const double& radius );
};

//------------------------------------------------------------------------------------------------------

struct Earth {
    // 6371229  -- IFS
    // 6367470  -- GRIB1
    // 6378137  -- WGS84 semi-major axis
    static constexpr double radiusInMeters() { return 6371229.; }
    static constexpr double radiusInKm() { return radiusInMeters() / 1e3; }

    static constexpr double areaInSqMeters() { return 4. * M_PI * radiusInMeters() * radiusInMeters(); }
    static constexpr double areaInSqKm() { return 4. * M_PI * radiusInKm() * radiusInKm(); }

    // Great-circle central angle between two points, in radians
    static double centralAngle( const PointLonLat&, const PointLonLat& );
    static double centralAngle( const PointXYZ&, const PointXYZ&, const double& radius = radiusInMeters() );

    // Great-circle distance between two points
    static double distanceInMeters( const PointLonLat&, const PointLonLat&, const double& radius = radiusInMeters() );
    static double distanceInMeters( const PointXYZ&, const PointXYZ&, const double& radius = radiusInMeters() );

    // Great-circle intermediate position provided two circle points and
    // longitude, in degrees
    static void greatCircleLatitudeGivenLongitude( const PointLonLat&, const PointLonLat&, PointLonLat& );

    // Convert geodetic coordinates to geocentric Cartesian (ECEF: Earth-centered,
    // Earth-fixed)
    static PointXYZ convertGeodeticToGeocentric( const PointLonLat&, const double& radius = radiusInMeters(),
                                                 const double& height = 0. );

    // Convert geocentric Cartesian (ECEF: Earth-centered, Earth-fixed) to
    // geodetic coordinates
    static PointLonLat convertGeocentricToGeodetic( const PointXYZ&, const double& radius = radiusInMeters() );
};

//------------------------------------------------------------------------------------------------------

}  // namespace util
}  // namespace atlas
