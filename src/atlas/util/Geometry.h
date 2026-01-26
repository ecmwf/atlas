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

#include <memory>
#include <sstream>
#include <string>
#include <utility>

#include "atlas/runtime/Exception.h"
#include "atlas/util/Earth.h"
#include "atlas/util/Object.h"
#include "atlas/util/ObjectHandle.h"
#include "atlas/util/Point.h"

namespace atlas {
namespace geometry {

//------------------------------------------------------------------------------------------------------

namespace detail {

// TODO: move greatCircleCourse to eckit::geometry::Sphere

/// @brief   Calculate great-cricle course between points
///
/// @details Calculates the direction (clockwise from north) of a great-circle
///          arc between lonLat1 and lonLat2. Returns the direction of the arc
///          at lonLat1 (first) and lonLat2 (second). Angle is normalised to the
///          range of atan2 (usually (-180, 180]). All input and output values
///          are in units of degrees.
/// @ref     https://en.wikipedia.org/wiki/Great-circle_navigation
std::pair<double, double> greatCircleCourse(const Point2& lonLat1, const Point2& lonLat2);

//------------------------------------------------------------------------------------------------------

class GeometryBase : public util::Object {
public:
    virtual ~GeometryBase()                                           = default;
    virtual void lonlat2xyz(const Point2&, Point3&) const             = 0;
    virtual void xyz2lonlat(const Point3&, Point2&) const             = 0;
    virtual double distance(const Point2& p1, const Point2& p2) const = 0;
    virtual double distance(const Point3& p1, const Point3& p2) const = 0;
    virtual double radius() const                                     = 0;
    virtual double area() const                                       = 0;
    virtual std::pair<double, double> greatCircleCourse(const Point2& p1, const Point2& p2) const = 0;

    Point3 xyz(const Point2& lonlat) const {
        Point3 xyz;
        lonlat2xyz(lonlat, xyz);
        return xyz;
    }
    Point2 lonlat(const Point3& xyz) const {
        Point2 lonlat;
        xyz2lonlat(xyz, lonlat);
        return lonlat;
    }
};

//------------------------------------------------------------------------------------------------------

template <typename SphereT>
class GeometrySphereT : public GeometryBase {
public:
    void lonlat2xyz(const Point2& lonlat, Point3& xyz) const override {
#if ATLAS_ECKIT_VERSION_AT_LEAST(1, 24, 0)
        SphereT::convertSphericalToCartesian(lonlat, xyz, 0., true);
#else
        SphereT::convertSphericalToCartesian(lonlat, xyz);
#endif
    }
    void xyz2lonlat(const Point3& xyz, Point2& lonlat) const override {
        SphereT::convertCartesianToSpherical(xyz, lonlat);
    }
    double distance(const Point2& p1, const Point2& p2) const override { return SphereT::distance(p1, p2); }
    double distance(const Point3& p1, const Point3& p2) const override { return SphereT::distance(p1, p2); }
    double radius() const override { return SphereT::radius(); }
    double area() const override { return SphereT::area(); }
    std::pair<double, double> greatCircleCourse(const Point2& p1, const Point2& p2) const override {
        return atlas::geometry::detail::greatCircleCourse(p1, p2);
    }
};

class GeometrySphere : public GeometryBase {
    using Sphere = eckit::geometry::Sphere;

public:
    GeometrySphere(double radius): radius_(radius) {}
    void lonlat2xyz(const Point2& lonlat, Point3& xyz) const override;
    void xyz2lonlat(const Point3& xyz, Point2& lonlat) const override;
    double distance(const Point2& p1, const Point2& p2) const override { return Sphere::distance(radius_, p1, p2); }
    double distance(const Point3& p1, const Point3& p2) const override { return Sphere::distance(radius_, p1, p2); }
    double radius() const override { return radius_; }
    double area() const override { return Sphere::area(radius_); }
    std::pair<double, double> greatCircleCourse(const Point2& p1, const Point2& p2) const override {
        return atlas::geometry::detail::greatCircleCourse(p1, p2);
    }

private:
    double radius_;
};

}  // namespace detail
}  // namespace geometry

//------------------------------------------------------------------------------------------------------

class Geometry : public util::ObjectHandle<geometry::detail::GeometryBase> {
public:
    using Handle::Handle;

    Geometry(): Handle(build<geometry::detail::GeometrySphereT<util::Earth>>()) {}
    Geometry(const std::string& name): Handle(build(name)) {}
    Geometry(const char* name): Handle(build(name)) {}
    Geometry(double radius): Handle(build<geometry::detail::GeometrySphere>(radius)) {}

    template <typename SphereT>
    Geometry(const SphereT&): Handle(build<SphereT>()) {}

    Point3 xyz(const Point2& lonlat) const { return get()->xyz(lonlat); }
    Point2 lonlat(const Point3& xyz) const { return get()->lonlat(xyz); }
    void xyz2lonlat(const Point3& xyz, Point2& lonlat) const { get()->xyz2lonlat(xyz, lonlat); }
    void lonlat2xyz(const Point2& lonlat, Point3& xyz) const { get()->lonlat2xyz(lonlat, xyz); }
    double distance(const Point2& p1, const Point2& p2) const { return get()->distance(p1, p2); }
    double distance(const Point3& p1, const Point3& p2) const { return get()->distance(p1, p2); }
    double radius() const { return get()->radius(); }
    double area() const { return get()->area(); }
    std::pair<double, double> greatCircleCourse(const Point2& p1, const Point2& p2) const { return get()->greatCircleCourse(p1, p2); }

protected:
    template <typename GeometryT, typename... Args>
    static Implementation* build(Args... args) {
        return new GeometryT(args...);
    }

    static Implementation* build(const std::string& name) {
        // Factory without self registration
        if (name == "Earth") {
            return build<geometry::detail::GeometrySphereT<util::Earth>>();
        }
        else if (name == "UnitSphere") {
            return build<geometry::detail::GeometrySphereT<eckit::geometry::UnitSphere>>();
        }
        else {
            ATLAS_THROW_EXCEPTION("name " << name << " is not a valid key for a Geometry");
        }
    }
};

//------------------------------------------------------------------------------------------------------

namespace geometry {
using Earth = Geometry;  // Sphere with util::Earth radius by default
class UnitSphere : public Geometry {
public:
    UnitSphere(): Geometry(/*radius*/ 1.) {}
};

}  // namespace geometry


// ------------------------------------------------------------------
// C wrapper interfaces to C++ routines

extern "C" {
Geometry::Implementation* atlas__Geometry__new_name(const char* name);
Geometry::Implementation* atlas__Geometry__new_radius(const double radius);
void atlas__Geometry__delete(Geometry::Implementation* This);
void atlas__Geometry__xyz2lonlat(Geometry::Implementation* This, const double x, const double y, const double z,
                                 double& lon, double& lat);
void atlas__Geometry__lonlat2xyz(Geometry::Implementation* This, const double lon, const double lat, double& x,
                                 double& y, double& z);
double atlas__Geometry__distance_lonlat(Geometry::Implementation* This, const double lon1, const double lat1,
                                        const double lon2, const double lat2);
double atlas__Geometry__distance_xyz(Geometry::Implementation* This, const double x1, const double y1, const double z1,
                                     const double x2, const double y2, const double z2);
double atlas__Geometry__radius(Geometry::Implementation* This);
double atlas__Geometry__area(Geometry::Implementation* This);
}

//------------------------------------------------------------------------------------------------------

}  // namespace atlas
