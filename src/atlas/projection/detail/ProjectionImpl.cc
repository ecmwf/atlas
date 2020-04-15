/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <algorithm>
#include <cstddef>
#include <memory>
#include <vector>

#include "eckit/types/FloatCompare.h"
#include "eckit/utils/Hash.h"

#include "ProjectionImpl.h"

#include "atlas/projection/detail/ProjectionFactory.h"
#include "atlas/runtime/Exception.h"
#include "atlas/util/Config.h"

namespace atlas {
namespace projection {
namespace detail {

// --------------------------------------------------------------------------------------------------------------------

namespace {

template <class T>
struct DerivateBuilder : public ProjectionImpl::DerivateFactory {
    using DerivateFactory::DerivateFactory;
    ProjectionImpl::Derivate* make( const ProjectionImpl& p, PointXY A, PointXY B, double h ) override {
        return new T( p, A, B, h );
    }
};

struct DerivateForwards final : ProjectionImpl::Derivate {
    using Derivate::Derivate;
    PointLonLat d( PointXY P ) const override {
        PointXY A( xy2lonlat( P ) );
        PointXY B( xy2lonlat( PointXY::add( P, H_ ) ) );
        return PointXY::div( PointXY::sub( B, A ), normH_ );
    }
};

struct DerivateBackwards final : ProjectionImpl::Derivate {
    using Derivate::Derivate;
    PointLonLat d( PointXY P ) const override {
        PointXY A( xy2lonlat( PointXY::sub( P, H_ ) ) );
        PointXY B( xy2lonlat( P ) );
        return PointXY::div( PointXY::sub( B, A ), normH_ );
    }
};

struct DerivateCentral final : ProjectionImpl::Derivate {
    DerivateCentral( const ProjectionImpl& p, PointXY A, PointXY B, double h ) :
        Derivate( p, A, B, h ), H2_{PointXY::mul( H_, 0.5 )} {}
    const PointXY H2_;
    PointLonLat d( PointXY P ) const override {
        PointXY A( xy2lonlat( PointXY::sub( P, H2_ ) ) );
        PointXY B( xy2lonlat( PointXY::add( P, H2_ ) ) );
        return PointXY::div( PointXY::sub( B, A ), normH_ );
    }
};

}  // namespace

ProjectionImpl::Derivate::Derivate( const ProjectionImpl& p, PointXY A, PointXY B, double h ) :
    projection_( p ),
    H_{PointXY::mul( PointXY::normalize( PointXY::sub( B, A ) ), h )},
    normH_( PointXY::norm( H_ ) ) {}

ProjectionImpl::Derivate::~Derivate() = default;

ProjectionImpl::Derivate* ProjectionImpl::DerivateFactory::build( const std::string& type, const ProjectionImpl& p,
                                                                  PointXY A, PointXY B, double h ) {
    ATLAS_ASSERT( 0. < h );

    // force_link
    static DerivateBuilder<DerivateForwards> __derivate1( "forwards" );
    static DerivateBuilder<DerivateBackwards> __derivate2( "backwards" );
    static DerivateBuilder<DerivateCentral> __derivate3( "central" );

    if ( A.distance2( B ) < h * h ) {
        struct DerivateDegenerate final : Derivate {
            using Derivate::Derivate;
            PointLonLat d( PointXY ) const override { return {}; }
        };
        return new DerivateDegenerate( p, A, B, h );
    }

    auto factory = get( type );
    return factory->make( p, A, B, h );
}

ProjectionImpl::DerivateFactory::~DerivateFactory() = default;

// --------------------------------------------------------------------------------------------------------------------

ProjectionImpl::BoundLonLat::operator RectangularLonLatDomain() const {
    return RectangularLonLatDomain( {min_[0], max_[0]}, {min_[1], max_[1]} );
}

bool ProjectionImpl::BoundLonLat::crossesDateLine( bool yes ) {
    if ( ( crossesDateLine_ = crossesDateLine_ || yes ) ) {
        max_.lon() = min_.lon() + 360.;
    }
    return crossesDateLine_;
}

bool ProjectionImpl::BoundLonLat::includesNorthPole( bool yes ) {
    if ( ( includesNorthPole_ = includesNorthPole_ || yes ) ) {
        max_.lat() = 90.;
    }
    crossesDateLine( includesNorthPole_ );
    return includesNorthPole_;
}

bool ProjectionImpl::BoundLonLat::includesSouthPole( bool yes ) {
    if ( ( includesSouthPole_ = includesSouthPole_ || yes ) ) {
        min_.lat() = -90.;
    }
    crossesDateLine( includesSouthPole_ );
    return includesSouthPole_;
}

void ProjectionImpl::BoundLonLat::extend( PointLonLat p, PointLonLat eps ) {
    ATLAS_ASSERT( 0. <= eps.lon() && 0. <= eps.lat() );

    auto sub = PointLonLat::sub( p, eps );
    auto add = PointLonLat::add( p, eps );
    min_     = first_ ? sub : PointLonLat::componentsMin( min_, sub );
    max_     = first_ ? add : PointLonLat::componentsMax( max_, add );
    first_   = false;

    min_.lat() = std::max( min_.lat(), -90. );
    max_.lat() = std::min( max_.lat(), 90. );
    max_.lon() = std::min( max_.lon(), min_.lon() + 360. );
    ATLAS_ASSERT( min_.lon() <= max_.lon() && min_.lat() <= max_.lat() );

    includesSouthPole( eckit::types::is_approximately_equal( min_.lat(), -90. ) );
    includesNorthPole( eckit::types::is_approximately_equal( max_.lat(), 90. ) );
    crossesDateLine( eckit::types::is_approximately_equal( max_.lon() - min_.lon(), 360. ) );
}

// --------------------------------------------------------------------------------------------------------------------

const ProjectionImpl* ProjectionImpl::create( const eckit::Parametrisation& p ) {
    std::string projectionType;
    if ( p.get( "type", projectionType ) ) {
        return ProjectionFactory::build( projectionType, p );
    }

    // should return error here
    throw_Exception( "type missing in Params", Here() );
}

PointXYZ ProjectionImpl::xyz( const PointLonLat& lonlat ) const {
    atlas::PointXYZ xyz;
    atlas::util::Earth::convertSphericalToCartesian( lonlat, xyz );
    return xyz;
}

RectangularLonLatDomain ProjectionImpl::lonlatBoundingBox( const Domain& domain ) const {
    using eckit::types::is_strictly_greater;


    // 0. setup

    if ( domain.global() ) {
        return domain;
    }
    RectangularDomain rect( domain );
    ATLAS_ASSERT( rect );

    constexpr double h     = 0.001;
    constexpr size_t Niter = 100;


    // 1. determine box from projected corners

    const std::vector<PointXY> corners{
        {rect.xmin(), rect.ymax()}, {rect.xmax(), rect.ymax()}, {rect.xmax(), rect.ymin()}, {rect.xmin(), rect.ymin()}};

    BoundLonLat bounds;
    for ( auto& p : corners ) {
        bounds.extend( lonlat( p ), PointLonLat{h, h} );
    }


    // 2. locate latitude extrema by checking if poles are included (in the un-projected frame) and if not, find extrema
    // not at the corners by refining iteratively

    for ( size_t i = 0; i < corners.size(); ++i ) {
        if ( !bounds.includesNorthPole() || !bounds.includesSouthPole() ) {
            PointXY A = corners[i];
            PointXY B = corners[( i + 1 ) % corners.size()];

            std::unique_ptr<Derivate> derivate( DerivateFactory::build( "central", *this, A, B ) );
            double dAdy = derivate->d( A ).lat();
            double dBdy = derivate->d( B ).lat();

            if ( !is_strictly_greater( dAdy * dBdy, 0. ) ) {
                for ( size_t cnt = 0; cnt < Niter; ++cnt ) {
                    PointXY M   = PointXY::middle( A, B );
                    double dMdy = derivate->d( M ).lat();
                    if ( is_strictly_greater( dAdy * dMdy, 0. ) ) {
                        A    = M;
                        dAdy = dMdy;
                    }
                    else if ( is_strictly_greater( dBdy * dMdy, 0. ) ) {
                        B    = M;
                        dBdy = dMdy;
                    }
                    else {
                        break;
                    }
                }

                // update extrema, extended by 'a small amount' (arbitrary)
                bounds.extend( lonlat( PointXY::middle( A, B ) ), PointLonLat{0, h} );
            }
        }
    }


    // 3. locate latitude extrema not at the corners by refining iteratively

    for ( size_t i = 0; i < corners.size(); ++i ) {
        if ( !bounds.crossesDateLine() ) {
            PointXY A = corners[i];
            PointXY B = corners[( i + 1 ) % corners.size()];

            std::unique_ptr<Derivate> derivate( DerivateFactory::build( "central", *this, A, B ) );
            double dAdx = derivate->d( A ).lon();
            double dBdx = derivate->d( B ).lon();

            if ( !is_strictly_greater( dAdx * dBdx, 0. ) ) {
                for ( size_t cnt = 0; cnt < Niter; ++cnt ) {
                    PointXY M   = PointXY::middle( A, B );
                    double dMdx = derivate->d( M ).lon();
                    if ( is_strictly_greater( dAdx * dMdx, 0. ) ) {
                        A    = M;
                        dAdx = dMdx;
                    }
                    else if ( is_strictly_greater( dBdx * dMdx, 0. ) ) {
                        B    = M;
                        dBdx = dMdx;
                    }
                    else {
                        break;
                    }
                }

                // update extrema, extended by 'a small amount' (arbitrary)
                bounds.extend( lonlat( PointXY::middle( A, B ) ), PointLonLat{h, 0} );
            }
        }
    }


    // 4. return bounding box
    return bounds;
}

//---------------------------------------------------------------------------------------------------------------------

Rotated::Rotated( const PointLonLat& south_pole, double rotation_angle ) :
    util::Rotation( south_pole, rotation_angle ) {}

Rotated::Rotated( const eckit::Parametrisation& p ) : util::Rotation( p ) {}

void Rotated::spec( Spec& s ) const {
    std::vector<double> npole{northPole().lon(), northPole().lat()};
    std::vector<double> spole{southPole().lon(), southPole().lat()};
    s.set( "north_pole", npole );
    s.set( "south_pole", spole );
    s.set( "rotation_angle", rotationAngle() );
}

void Rotated::hash( eckit::Hash& hsh ) const {
    hsh.add( "rotated" );
    hsh.add( southPole().lon() );
    hsh.add( southPole().lat() );
    hsh.add( rotationAngle() );
}

}  // namespace detail
}  // namespace projection
}  // namespace atlas
