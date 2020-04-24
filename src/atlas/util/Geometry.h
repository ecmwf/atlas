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
#include <string>

#include "atlas/runtime/Exception.h"
#include "atlas/util/Earth.h"
#include "atlas/util/Point.h"

namespace atlas {
namespace util {

struct GeometryBase {
    virtual ~GeometryBase()                                 = default;
    virtual void lonlat2xyz( const Point2&, Point3& ) const = 0;
    virtual void xyz2lonlat( const Point3&, Point2& ) const = 0;

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

template <typename SphereT>
struct GeometrySphereT : public GeometryBase {
    void lonlat2xyz( const Point2& lonlat, Point3& xyz ) const override {
        SphereT::convertSphericalToCartesian( lonlat, xyz );
    }
    void xyz2lonlat( const Point3& xyz, Point2& lonlat ) const override {
        SphereT::convertCartesianToSpherical( xyz, lonlat );
    }
};

struct GeometrySphere : public GeometryBase {
    GeometrySphere( double radius ) : radius_( radius ) {}
    void lonlat2xyz( const Point2& lonlat, Point3& xyz ) const override {
        eckit::geometry::Sphere::convertSphericalToCartesian( radius_, lonlat, xyz );
    }
    void xyz2lonlat( const Point3& xyz, Point2& lonlat ) const override {
        eckit::geometry::Sphere::convertCartesianToSpherical( radius_, xyz, lonlat );
    }
    double radius_;
};

struct Geometry {
    Geometry() { build(); }

    Geometry( const std::string& name ) { build( name ); }

    Geometry( const std::shared_ptr<GeometryBase>& conversion ) :
        shared_impl_( conversion ), impl_( shared_impl_.get() ) {}

    Geometry( const Geometry& other ) : shared_impl_( other.shared_impl_ ), impl_( shared_impl_.get() ) {}

    Point3 xyz( const Point2& lonlat ) const { return impl_->xyz( lonlat ); }
    Point2 lonlat( const Point3& xyz ) const { return impl_->lonlat( xyz ); }
    void xyz2lonlat( const Point3& xyz, Point2& lonlat ) const { return impl_->xyz2lonlat( xyz, lonlat ); }
    void lonlat2xyz( const Point2& lonlat, Point3& xyz ) const { return impl_->lonlat2xyz( lonlat, xyz ); }

private:
    std::shared_ptr<GeometryBase> shared_impl_;
    GeometryBase* impl_;

    template <typename GeometryT, typename... Args>
    void build( Args... args ) {
        shared_impl_ = std::make_shared<GeometryT>( args... );
        impl_        = shared_impl_.get();
    }

    void build( const std::string& name = "Earth" ) {
        // Factory without self registration
        if ( name == "Earth" ) {
            build<GeometrySphereT<Earth>>();
        }
        else if ( name == "UnitSphere" ) {
            build<GeometrySphereT<eckit::geometry::UnitSphere>>();
        }
        else {
            ATLAS_THROW_EXCEPTION( "name " + name + " is not a valid key for a Geometry" );
        }
    }
};

struct Sphere : public Geometry {
    using Geometry::Geometry;
    Sphere( double radius = Earth::radius() ) : Geometry( std::make_shared<GeometrySphere>( radius ) ) {}
};

//------------------------------------------------------------------------------------------------------

}  // namespace util
}  // namespace atlas
