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

#include "eckit/geometry/Point2.h"
#include "eckit/geometry/Point3.h"
#include "eckit/geometry/UnitSphere.h"

//------------------------------------------------------------------------------------------------------

namespace atlas {
namespace util {

//------------------------------------------------------------------------------------------------------

using eckit::geometry::Point2;
using eckit::geometry::Point3;
using eckit::geometry::UnitSphere;

// ------------------------------------------------------------------
// C wrapper interfaces to C++ routines

extern "C" {
UnitSphere* atlas__UnitSphere__new();
void atlas__UnitSphere__delete( UnitSphere* This );
double atlas__UnitSphere__radius( UnitSphere* This );
double atlas__UnitSphere__central_angle_2( UnitSphere* This, const Point2* Alonlat, const Point2* Blonlat );
double atlas__UnitSphere__central_angle_3( UnitSphere* This, const Point3* A, const Point3* B );
double atlas__UnitSphere__distance_2( UnitSphere* This, const Point2* Alonlat, const Point2* Blonlat );
double atlas__UnitSphere__distance_3( UnitSphere* This, const Point3* A, const Point3* B );
double atlas__UnitSphere__area( UnitSphere* This );
double atlas__UnitSphere__area_wn_es( UnitSphere* This, const Point2* WestNorth, const Point2* EastSouth );
double atlas__UnitSphere__great_circle_latitude_given_longitude( UnitSphere* This, const Point2* Alonlat, const Point2* Blonlat, const double& Clon );
void atlas__UnitSphere__great_circle_longitude_given_latitude( UnitSphere* This, const Point2* Alonlat, const Point2* Blonlat, const double& Clat, double& Clon1, double& Clon2 );
void atlas__UnitSphere__convert_spherical_to_cartesian( UnitSphere* This, const Point2* Alonlat, Point3* B );
void atlas__UnitSphere__convert_cartesian_to_spherical( UnitSphere* This, const Point3* A, Point2* Blonlat ); 
}

//------------------------------------------------------------------------------------------------------

}  // namespace util
}  // namespace atlas
