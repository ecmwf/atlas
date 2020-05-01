/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immGeometryies
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "eckit/geometry/Point2.h"
#include "eckit/geometry/Point3.h"

#include "atlas/util/Geometry.h"
#include "atlas/runtime/Exception.h"

namespace atlas {
namespace geometry {

// ------------------------------------------------------------------
// C wrapper interfaces to C++ routines
Geometry* atlas__Geometry__new_name( const char* name ) {
    return new Geometry(std::string(name));
}
Geometry* atlas__Geometry__new_radius( const double radius ) {
    return new Geometry(radius);
}
void atlas__Geometry__delete( Geometry* This ) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialised atlas_Geometry" );
    delete This;
}
void atlas__Geometry__xyz2lonlat( Geometry* This, const Point3* xyz, Point2* lonlat ) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialised atlas_Geometry" );
    return This->xyz2lonlat(*xyz, *lonlat);
}
void atlas__Geometry__lonlat2xyz( Geometry* This, const Point2* lonlat, Point3* xyz ) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialised atlas_Geometry" );
    return This->lonlat2xyz(*lonlat, *xyz);
}
double atlas__Geometry__distance_2( Geometry* This, const Point2* p1, const Point2* p2 ) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialised atlas_Geometry" );
    return This->distance(*p1, *p2);
}
double atlas__Geometry__distance_3( Geometry* This, const Point3* p1, const Point3* p2 ) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialised atlas_Geometry" );
    return This->distance(*p1, *p2);
}
double atlas__Geometry__radius( Geometry* This) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialised atlas_Geometry" );
    return This->radius();
}
double atlas__Geometry__area( Geometry* This) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialised atlas_Geometry" );
    return This->area();
}

// ------------------------------------------------------------------

}  // namespace util
}  // namespace atlas
