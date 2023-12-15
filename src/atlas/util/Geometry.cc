/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immGeometryies
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/util/Geometry.h"

#include <cmath>

#include "eckit/geometry/Point2.h"
#include "eckit/geometry/Point3.h"

#include "atlas/library/config.h"
#include "atlas/runtime/Exception.h"
#include "atlas/util/Constants.h"

namespace atlas {

namespace geometry {
namespace detail {
void GeometrySphere::lonlat2xyz(const Point2& lonlat, Point3& xyz) const {
#if ATLAS_ECKIT_VERSION_AT_LEAST(1, 24, 0)
    Sphere::convertSphericalToCartesian(radius_, lonlat, xyz, 0., true);
#else
    Sphere::convertSphericalToCartesian(radius_, lonlat, xyz);
#endif
}
void GeometrySphere::xyz2lonlat(const Point3& xyz, Point2& lonlat) const {
    Sphere::convertCartesianToSpherical(radius_, xyz, lonlat);
}


/// @brief   Calculate great-cricle course between points
///
/// @details Calculates the direction (clockwise from north) of a great-circle
///          arc between lonLat1 and lonLat2. Returns the direction of the arc
///          at lonLat1 (first) and lonLat2 (second). Angle is normalised to the
///          range of atan2 (usually (-180, 180]). All input and output values
///          are in units of degrees.
/// @ref     https://en.wikipedia.org/wiki/Great-circle_navigation
///
std::pair<double, double> greatCircleCourse(const Point2& lonLat1,
                                                   const Point2& lonLat2) {

  const auto lambda1 = lonLat1[0] * util::Constants::degreesToRadians();
  const auto lambda2 = lonLat2[0] * util::Constants::degreesToRadians();
  const auto phi1 = lonLat1[1] * util::Constants::degreesToRadians();
  const auto phi2 = lonLat2[1] * util::Constants::degreesToRadians();

  const auto sinLambda12 = std::sin(lambda2 - lambda1);
  const auto cosLambda12 = std::cos(lambda2 - lambda1);
  const auto sinPhi1 = std::sin(phi1);
  const auto sinPhi2 = std::sin(phi2);
  const auto cosPhi1 = std::cos(phi1);
  const auto cosPhi2 = std::cos(phi2);

  const auto alpha1 =
      std::atan2(cosPhi2 * sinLambda12,
                 cosPhi1 * sinPhi2 - sinPhi1 * cosPhi2 * cosLambda12);

  const auto alpha2 =
      std::atan2(cosPhi1 * sinLambda12,
                 -cosPhi2 * sinPhi1 + sinPhi2 * cosPhi1 * cosLambda12);

  return std::make_pair(alpha1 * util::Constants::radiansToDegrees(),
                        alpha2 * util::Constants::radiansToDegrees());
};


} // namespace detail
} // namespace geometry

extern "C" {
// ------------------------------------------------------------------
// C wrapper interfaces to C++ routines
Geometry::Implementation* atlas__Geometry__new_name(const char* name) {
    Geometry::Implementation* geometry;
    {
        Geometry handle{std::string{name}};
        geometry = handle.get();
        geometry->attach();
    }
    geometry->detach();
    return geometry;
}
geometry::detail::GeometryBase* atlas__Geometry__new_radius(const double radius) {
    Geometry::Implementation* geometry;
    {
        Geometry handle{radius};
        geometry = handle.get();
        geometry->attach();
    }
    geometry->detach();
    return geometry;
}
void atlas__Geometry__delete(Geometry::Implementation* This) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_Geometry");
    delete This;
}
void atlas__Geometry__xyz2lonlat(Geometry::Implementation* This, const double x, const double y, const double z,
                                 double& lon, double& lat) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_Geometry");
    PointLonLat lonlat;
    This->xyz2lonlat(PointXYZ{x, y, z}, lonlat);
    lon = lonlat.lon();
    lat = lonlat.lat();
}
void atlas__Geometry__lonlat2xyz(Geometry::Implementation* This, const double lon, const double lat, double& x,
                                 double& y, double& z) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_Geometry");
    PointXYZ xyz;
    This->lonlat2xyz(PointLonLat{lon, lat}, xyz);
    x = xyz.x();
    y = xyz.y();
    z = xyz.z();
}
double atlas__Geometry__distance_lonlat(Geometry::Implementation* This, const double lon1, const double lat1,
                                        const double lon2, const double lat2) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_Geometry");
    return This->distance(PointLonLat{lon1, lat1}, PointLonLat{lon2, lat2});
}
double atlas__Geometry__distance_xyz(Geometry::Implementation* This, const double x1, const double y1, const double z1,
                                     const double x2, const double y2, const double z2) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_Geometry");
    return This->distance(PointXYZ{x1, y1, z1}, PointXYZ{x2, y2, z2});
}
double atlas__Geometry__radius(Geometry::Implementation* This) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_Geometry");
    return This->radius();
}
double atlas__Geometry__area(Geometry::Implementation* This) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_Geometry");
    return This->area();
}
}
// ------------------------------------------------------------------

}  // namespace atlas
