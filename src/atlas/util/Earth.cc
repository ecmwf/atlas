/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "eckit/geometry/Point2.h"
#include "eckit/geometry/Point3.h"

#include "atlas/util/Earth.h"
#include "atlas/runtime/Exception.h"

namespace atlas {
namespace util {

// ------------------------------------------------------------------
// C wrapper interfaces to C++ routines

Earth* atlas__Earth__new() {
    return new Earth();
}
void atlas__Earth__delete( Earth* This ) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialised atlas_Earth" );
    delete This;
}
double atlas__Earth__radius( Earth* This) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialised atlas_Earth" );
    return This->radius();
}
double atlas__Earth__central_angle_2( Earth* This, const Point2* Alonlat, const Point2* Blonlat ) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialised atlas_Earth" );
    return This->centralAngle(*Alonlat, *Blonlat);
}
double atlas__Earth__central_angle_3( Earth* This, const Point3* A, const Point3* B ) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialised atlas_Earth" );
    return This->centralAngle(*A, *B);
}
double atlas__Earth__distance_2( Earth* This, const Point2* Alonlat, const Point2* Blonlat ) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialised atlas_Earth" );
    return This->distance(*Alonlat, *Blonlat);
}
double atlas__Earth__distance_3( Earth* This, const Point3* A, const Point3* B ) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialised atlas_Earth" );
    return This->distance(*A, *B);
}
double atlas__Earth__area( Earth* This) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialised atlas_Earth" );
    return This->area();
}
double atlas__Earth__area_wn_es( Earth* This, const Point2* WestNorth, const Point2* EastSouth ) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialised atlas_Earth" );
    return This->area(*WestNorth, *EastSouth);
}
double atlas__Earth__great_circle_latitude_given_longitude( Earth* This, const Point2* Alonlat, const Point2* Blonlat, const double& Clon ) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialised atlas_Earth" );
    return This->greatCircleLatitudeGivenLongitude(*Alonlat, *Blonlat, Clon);
}
void atlas__Earth__great_circle_longitude_given_latitude( Earth* This, const Point2* Alonlat, const Point2* Blonlat, const double& Clat, double& Clon1, double& Clon2 ) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialised atlas_Earth" );
    return This->greatCircleLongitudeGivenLatitude(*Alonlat, *Blonlat, Clat, Clon1, Clon2);
}
void atlas__Earth__convert_spherical_to_cartesian( Earth* This, const Point2* Alonlat, Point3* B ) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialised atlas_Earth" );
    return This->convertSphericalToCartesian(*Alonlat, *B);
}
void atlas__Earth__convert_cartesian_to_spherical( Earth* This, const Point3* A, Point2* Blonlat ) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialised atlas_Earth" );
    return This->convertCartesianToSpherical(*A, *Blonlat);
}

// ------------------------------------------------------------------

}  // namespace util
}  // namespace atlas
