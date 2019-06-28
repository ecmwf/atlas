/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <memory>

#include "eckit/thread/AutoLock.h"
#include "eckit/types/FloatCompare.h"
#include "eckit/utils/Hash.h"

#include "atlas/interpolation/element/Quad3D.h"
#include "atlas/interpolation/method/Ray.h"
#include "atlas/projection/detail/LonLatProjection.h"
#include "atlas/projection/detail/ProjectionFactory.h"
#include "atlas/projection/detail/ProjectionImpl.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"
#include "atlas/util/Config.h"
#include "atlas/util/Earth.h"

namespace atlas {
namespace projection {
namespace detail {

// --------------------------------------------------------------------------------------------------------------------

namespace {

template <class T>
struct DerivateBuilder : public ProjectionImpl::DerivateFactory {
    DerivateBuilder( const std::string& type ) : DerivateFactory( type ) {}
    ProjectionImpl::Derivate* make( const ProjectionImpl& p, PointXY A, PointXY B, double h ) {
        return new T( p, A, B, h );
    }
};

static pthread_once_t once = PTHREAD_ONCE_INIT;
static eckit::Mutex* mtx   = nullptr;

static std::map<std::string, ProjectionImpl::DerivateFactory*>* m = nullptr;
static void init() {
    mtx = new eckit::Mutex();
    m   = new std::map<std::string, ProjectionImpl::DerivateFactory*>();
}

}  // (anonymous namespace)

ProjectionImpl::Derivate* ProjectionImpl::DerivateFactory::build( const std::string& type, const ProjectionImpl& p,
                                                                  PointXY A, PointXY B, double h ) {
    pthread_once( &once, init );
    eckit::AutoLock<eckit::Mutex> __lock( mtx );

    ATLAS_ASSERT( 0. < h );
    if ( A.distance2( B ) < h * h ) {
        struct DerivateDegenerate final : Derivate {
            using Derivate::Derivate;
            PointLonLat d( PointXY ) const override { return {}; }
        };
        return new DerivateDegenerate( p, A, B, h );
    }

    auto j = m->find( type );
    if ( j == m->end() ) {
        list( Log::error() << "DerivateFactory: unknown '" << type << "', choices are: " );
        throw_Exception( "DerivateFactory: unknown '" + type + "'", Here() );
    }
    return ( *j ).second->make( p, A, B, h );
}

void ProjectionImpl::DerivateFactory::list( std::ostream& out ) {
    pthread_once( &once, init );
    eckit::AutoLock<eckit::Mutex> __lock( mtx );

    const char* sep = "";
    for ( const auto& j : *m ) {
        out << sep << j.first;
        sep = ", ";
    }
}

ProjectionImpl::DerivateFactory::DerivateFactory( const std::string& type ) {
    pthread_once( &once, init );
    eckit::AutoLock<eckit::Mutex> __lock( mtx );

    if ( m->find( type ) != m->end() ) {
        throw_Exception( "DerivateFactory: duplicate '" + type + "'", Here() );
    }
    ( *m )[type] = this;
}

ProjectionImpl::DerivateFactory::~DerivateFactory() = default;

// --------------------------------------------------------------------------------------------------------------------

ProjectionImpl::BoundLonLat::operator Domain() const {
    return RectangularDomain( {min_[0], max_[0]}, {min_[1], max_[1]}, "degrees" );
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
    ATLAS_ASSERT( PointLonLat{} < eps );
    if ( first_ ) {
        min_   = PointLonLat::sub( p, eps );
        max_   = PointLonLat::add( p, eps );
        first_ = false;
        return;
    }

    min_ = PointLonLat::componentsMin( min_, PointLonLat::sub( p, eps ) );
    max_ = PointLonLat::componentsMax( max_, PointLonLat::add( p, eps ) );

    min_.lat() = std::max( min_.lat(), -90. );
    max_.lat() = std::min( max_.lat(), 90. );
    max_.lon() = std::min( max_.lon(), min_.lon() + 360. );
    ATLAS_ASSERT( min_ < max_ );

    includesSouthPole( eckit::types::is_approximately_equal( min_.lat(), -90. ) );
    includesNorthPole( eckit::types::is_approximately_equal( max_.lat(), 90. ) );
    crossesDateLine( eckit::types::is_approximately_equal( max_.lon() - min_.lon(), 360. ) );
}

// --------------------------------------------------------------------------------------------------------------------

ProjectionImpl::Derivate::Derivate( const ProjectionImpl& p, PointXY A, PointXY B, double h ) :
    projection_( p ),
    H_{PointXY::mul( PointXY::normalize( PointXY::sub( B, A ) ), h )},
    normH_( PointXY::norm( H_ ) ) {}

ProjectionImpl::Derivate::~Derivate() = default;

namespace {

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
        Derivate( p, A, B, h ),
        H2_{PointXY::mul( H_, 0.5 )} {}
    const PointXY H2_;
    PointLonLat d( PointXY P ) const override {
        PointXY A( xy2lonlat( PointXY::sub( P, H2_ ) ) );
        PointXY B( xy2lonlat( PointXY::add( P, H2_ ) ) );
        return PointXY::div( PointXY::sub( B, A ), normH_ );
    }
};

DerivateBuilder<DerivateForwards> __derivate1( "forwards" );
DerivateBuilder<DerivateBackwards> __derivate2( "backwards" );
DerivateBuilder<DerivateCentral> __derivate3( "central" );

}  // (anonymous namespace)

// --------------------------------------------------------------------------------------------------------------------

const ProjectionImpl* ProjectionImpl::create( const eckit::Parametrisation& p ) {
    std::string projectionType;
    if ( p.get( "type", projectionType ) ) {
        return ProjectionFactory::build( projectionType, p );
    }

    // should return error here
    throw_Exception( "type missing in Params", Here() );
}

Domain ProjectionImpl::boundingBox( const Domain& domain ) const {
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

    {
        auto xyz = [=]( const PointXY& p ) -> PointXYZ {
            PointXYZ r;
            util::Earth::convertSphericalToCartesian( this->lonlat( p ), r );
            return r;
        };

        using interpolation::element::Quad3D;
        using interpolation::method::Ray;
        const Quad3D quad( xyz( corners[0] ), xyz( corners[1] ), xyz( corners[2] ), xyz( corners[3] ) );

        PointXY NP{xy( {0., 90.} )};
        PointXY SP{xy( {0., -90.} )};

        bounds.includesNorthPole( quad.intersects( Ray( xyz( NP ) ) ) );
        bounds.includesNorthPole( quad.intersects( Ray( xyz( SP ) ) ) );
    }

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
