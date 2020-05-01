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

#include "atlas/runtime/Exception.h"
#include "atlas/util/Earth.h"
#include "atlas/util/Object.h"
#include "atlas/util/ObjectHandle.h"
#include "atlas/util/Point.h"

namespace atlas {
namespace geometry {

//------------------------------------------------------------------------------------------------------

namespace detail {

class GeometryBase : util::Object {
public:
    virtual ~GeometryBase()                                             = default;
    virtual void lonlat2xyz( const Point2&, Point3& ) const             = 0;
    virtual void xyz2lonlat( const Point3&, Point2& ) const             = 0;
    virtual double distance( const Point2& p1, const Point2& p2 ) const = 0;
    virtual double distance( const Point3& p1, const Point3& p2 ) const = 0;
    virtual double radius() const                                       = 0;
    virtual double area() const                                         = 0;

    Point3 xyz( const Point2& lonlat ) const {
        Point3 xyz;
        lonlat2xyz( lonlat, xyz );
        return xyz;
    }
    Point2 lonlat( const Point3& xyz ) const {
        Point2 lonlat;
        xyz2lonlat( xyz, lonlat );
        return lonlat;
    }
};

//------------------------------------------------------------------------------------------------------

template <typename SphereT>
class GeometrySphereT : public GeometryBase {
public:
    void lonlat2xyz( const Point2& lonlat, Point3& xyz ) const override {
        SphereT::convertSphericalToCartesian( lonlat, xyz );
    }
    void xyz2lonlat( const Point3& xyz, Point2& lonlat ) const override {
        SphereT::convertCartesianToSpherical( xyz, lonlat );
    }
    double distance( const Point2& p1, const Point2& p2 ) const override { return SphereT::distance( p1, p2 ); }
    double distance( const Point3& p1, const Point3& p2 ) const override { return SphereT::distance( p1, p2 ); }
    double radius() const override { return SphereT::radius(); }
    double area() const override { return SphereT::area(); }
};

class GeometrySphere : public GeometryBase {
    using Sphere = eckit::geometry::Sphere;

public:
    GeometrySphere( double radius ) : radius_( radius ) {}
    void lonlat2xyz( const Point2& lonlat, Point3& xyz ) const override {
        Sphere::convertSphericalToCartesian( radius_, lonlat, xyz );
    }
    void xyz2lonlat( const Point3& xyz, Point2& lonlat ) const override {
        Sphere::convertCartesianToSpherical( radius_, xyz, lonlat );
    }
    double distance( const Point2& p1, const Point2& p2 ) const override { return Sphere::distance( radius_, p1, p2 ); }
    double distance( const Point3& p1, const Point3& p2 ) const override { return Sphere::distance( radius_, p1, p2 ); }
    double radius() const override { return radius_; }
    double area() const override { return Sphere::area( radius_ ); }

private:
    double radius_;
};

}  // namespace detail
}  // namespace geometry

//------------------------------------------------------------------------------------------------------

class Geometry : util::ObjectHandle<geometry::detail::GeometryBase> {
public:
    using Handle::Handle;

    Geometry() : Handle( build<geometry::detail::GeometrySphereT<util::Earth>>() ) {}
    Geometry( const std::string& name ) : Handle( build( name ) ) {}
    Geometry( double radius ) : Handle( build<geometry::detail::GeometrySphere>( radius ) ) {}

    template <typename SphereT>
    Geometry( const SphereT& ) : Handle( build<SphereT>() ) {}

    Point3 xyz( const Point2& lonlat ) const { return get()->xyz( lonlat ); }
    Point2 lonlat( const Point3& xyz ) const { return get()->lonlat( xyz ); }
    void xyz2lonlat( const Point3& xyz, Point2& lonlat ) const { get()->xyz2lonlat( xyz, lonlat ); }
    void lonlat2xyz( const Point2& lonlat, Point3& xyz ) const { get()->lonlat2xyz( lonlat, xyz ); }
    double distance( const Point2& p1, const Point2& p2 ) const { return get()->distance( p1, p2 ); }
    double distance( const Point3& p1, const Point3& p2 ) const { return get()->distance( p1, p2 ); }
    double radius() const { return get()->radius(); }
    double area() const { return get()->area(); }

protected:
    template <typename GeometryT, typename... Args>
    static Implementation* build( Args... args ) {
        return new GeometryT( args... );
    }

    static Implementation* build( const std::string& name ) {
        // Factory without self registration
        if ( name == "Earth" ) {
            return build<geometry::detail::GeometrySphereT<util::Earth>>();
        }
        else if ( name == "UnitSphere" ) {
            return build<geometry::detail::GeometrySphereT<eckit::geometry::UnitSphere>>();
        }
        else {
            ATLAS_THROW_EXCEPTION( "name " << name << " is not a valid key for a Geometry" );
        }
    }
};

//------------------------------------------------------------------------------------------------------

namespace geometry {
using Earth      = Geometry;  // Sphere with util::Earth radius by default
using UnitSphere = Geometry( eckit::geometry::UnitSphere() );
}  // namespace geometry

//------------------------------------------------------------------------------------------------------

}  // namespace atlas
