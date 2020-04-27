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

#include "atlas/util/UnitSphere.h"
#include "atlas/runtime/Exception.h"

namespace atlas {
namespace util {

// ------------------------------------------------------------------
// C wrapper interfaces to C++ routines

UnitSphere* atlas__UnitSphere__new() {
    return new UnitSphere();
}
void atlas__UnitSphere__delete( UnitSphere* This ) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialised atlas_UnitSphere" );
    delete This;
}
double atlas__UnitSphere__radius( UnitSphere* This) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialised atlas_UnitSphere" );
    return This->radius();
}
double atlas__UnitSphere__central_angle_2( UnitSphere* This, const Point2* Alonlat, const Point2* Blonlat ) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialised atlas_UnitSphere" );
    return This->centralAngle(*Alonlat, *Blonlat);
}
double atlas__UnitSphere__central_angle_3( UnitSphere* This, const Point3* A, const Point3* B ) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialised atlas_UnitSphere" );
    return This->centralAngle(*A, *B);
}
double atlas__UnitSphere__distance_2( UnitSphere* This, const Point2* Alonlat, const Point2* Blonlat ) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialised atlas_UnitSphere" );
    return This->distance(*Alonlat, *Blonlat);
}
double atlas__UnitSphere__distance_3( UnitSphere* This, const Point3* A, const Point3* B ) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialised atlas_UnitSphere" );
    return This->distance(*A, *B);
}
double atlas__UnitSphere__area( UnitSphere* This) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialised atlas_UnitSphere" );
    return This->area();
}
double atlas__UnitSphere__area_wn_es( UnitSphere* This, const Point2* WestNorth, const Point2* EastSouth ) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialised atlas_UnitSphere" );
    return This->area(*WestNorth, *EastSouth);
}
double atlas__UnitSphere__great_circle_latitude_given_longitude( UnitSphere* This, const Point2* Alonlat, const Point2* Blonlat, const double& Clon ) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialised atlas_UnitSphere" );
    return This->greatCircleLatitudeGivenLongitude(*Alonlat, *Blonlat, Clon);
}
void atlas__UnitSphere__great_circle_longitude_given_latitude( UnitSphere* This, const Point2* Alonlat, const Point2* Blonlat, const double& Clat, double& Clon1, double& Clon2 ) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialised atlas_UnitSphere" );
    return This->greatCircleLongitudeGivenLatitude(*Alonlat, *Blonlat, Clat, Clon1, Clon2);
}
void atlas__UnitSphere__convert_spherical_to_cartesian( UnitSphere* This, const Point2* Alonlat, Point3* B ) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialised atlas_UnitSphere" );
    return This->convertSphericalToCartesian(*Alonlat, *B);
}
void atlas__UnitSphere__convert_cartesian_to_spherical( UnitSphere* This, const Point3* A, Point2* Blonlat ) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialised atlas_UnitSphere" );
    return This->convertCartesianToSpherical(*A, *Blonlat);
}

// ------------------------------------------------------------------

}  // namespace util
}  // namespace atlas
