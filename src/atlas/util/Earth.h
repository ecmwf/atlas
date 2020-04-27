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
#include "eckit/geometry/SphereT.h"
#include "eckit/geometry/UnitSphere.h"

//------------------------------------------------------------------------------------------------------

namespace atlas {
namespace util {

//------------------------------------------------------------------------------------------------------

using eckit::geometry::Point2;
using eckit::geometry::Point3;

//------------------------------------------------------------------------------------------------------

struct DatumIFS {
    static constexpr double radius() { return 6371229.; }
};

struct DatumGRIB1 {
    static constexpr double radius() { return 6367470.; }
};

struct DatumWGS84SemiMajorAxis {
    static constexpr double radius() { return 6378137.; }
};

//------------------------------------------------------------------------------------------------------

typedef eckit::geometry::SphereT<DatumIFS> Earth;

// ------------------------------------------------------------------
// C wrapper interfaces to C++ routines

extern "C" {
Earth* atlas__Earth__new();
void atlas__Earth__delete( Earth* This );
double atlas__Earth__radius( Earth* This );
double atlas__Earth__central_angle_2( Earth* This, const Point2* Alonlat, const Point2* Blonlat );
double atlas__Earth__central_angle_3( Earth* This, const Point3* A, const Point3* B );
double atlas__Earth__distance_2( Earth* This, const Point2* Alonlat, const Point2* Blonlat );
double atlas__Earth__distance_3( Earth* This, const Point3* A, const Point3* B );
double atlas__Earth__area( Earth* This );
double atlas__Earth__area_wn_es( Earth* This, const Point2* WestNorth, const Point2* EastSouth );
double atlas__Earth__great_circle_latitude_given_longitude( Earth* This, const Point2* Alonlat, const Point2* Blonlat, const double& Clon );
void atlas__Earth__great_circle_longitude_given_latitude( Earth* This, const Point2* Alonlat, const Point2* Blonlat, const double& Clat, double& Clon1, double& Clon2 );
void atlas__Earth__convert_spherical_to_cartesian( Earth* This, const Point2* Alonlat, Point3* B );
void atlas__Earth__convert_cartesian_to_spherical( Earth* This, const Point3* A, Point2* Blonlat ); 
}

//------------------------------------------------------------------------------------------------------

}  // namespace util
}  // namespace atlas
