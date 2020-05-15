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

#include "atlas/runtime/Exception.h"
#include "atlas/util/Geometry.h"

namespace atlas {

extern "C" {
// ------------------------------------------------------------------
// C wrapper interfaces to C++ routines
Geometry::Implementation* atlas__Geometry__new_name( const char* name ) {
    Geometry::Implementation* geometry;
    {
        Geometry handle{std::string{name}};
        geometry = handle.get();
        geometry->attach();
    }
    geometry->detach();
    return geometry;
}
geometry::detail::GeometryBase* atlas__Geometry__new_radius( const double radius ) {
    Geometry::Implementation* geometry;
    {
        Geometry handle{radius};
        geometry = handle.get();
        geometry->attach();
    }
    geometry->detach();
    return geometry;
}
void atlas__Geometry__delete( Geometry::Implementation* This ) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialised atlas_Geometry" );
    delete This;
}
void atlas__Geometry__xyz2lonlat( Geometry::Implementation* This, const double x, const double y, const double z,
                                  double& lon, double& lat ) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialised atlas_Geometry" );
    PointLonLat lonlat;
    This->xyz2lonlat( PointXYZ{x, y, z}, lonlat );
    lon = lonlat.lon();
    lat = lonlat.lat();
}
void atlas__Geometry__lonlat2xyz( Geometry::Implementation* This, const double lon, const double lat, double& x,
                                  double& y, double& z ) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialised atlas_Geometry" );
    PointXYZ xyz;
    This->lonlat2xyz( PointLonLat{lon, lat}, xyz );
    x = xyz.x();
    y = xyz.y();
    z = xyz.z();
}
double atlas__Geometry__distance_lonlat( Geometry::Implementation* This, const double lon1, const double lat1,
                                         const double lon2, const double lat2 ) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialised atlas_Geometry" );
    return This->distance( PointLonLat{lon1, lat1}, PointLonLat{lon2, lat2} );
}
double atlas__Geometry__distance_xyz( Geometry::Implementation* This, const double x1, const double y1, const double z1,
                                      const double x2, const double y2, const double z2 ) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialised atlas_Geometry" );
    return This->distance( PointXYZ{x1, y1, z1}, PointXYZ{x2, y2, z2} );
}
double atlas__Geometry__radius( Geometry::Implementation* This ) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialised atlas_Geometry" );
    return This->radius();
}
double atlas__Geometry__area( Geometry::Implementation* This ) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialised atlas_Geometry" );
    return This->area();
}
}
// ------------------------------------------------------------------

}  // namespace atlas
