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
    ProjectionImpl::Derivate* make( const ProjectionImpl& p, PointXY A, PointXY B, double h ) { return new T( p, A, B, h ); }
};

static pthread_once_t once = PTHREAD_ONCE_INIT;
static eckit::Mutex* mtx  = nullptr;

static std::map<std::string, ProjectionImpl::DerivateFactory*>* m = nullptr;
static void init() {
    mtx = new eckit::Mutex();
    m    = new std::map<std::string, ProjectionImpl::DerivateFactory*>();
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
    if ( p.get( "type", projectionType ) ) { return ProjectionFactory::build( projectionType, p ); }

    // should return error here
    throw_Exception( "type missing in Params", Here() );
}

Domain ProjectionImpl::boundingBox( const Domain& domain ) const {
    using eckit::types::is_approximately_lesser_or_equal;
    using eckit::types::is_strictly_greater;

    // 0. setup

    if ( domain.global() ) {
        return domain;
    }
    RectangularDomain rect( domain );
    ATLAS_ASSERT( rect );

    PointLonLat min;
    PointLonLat max;
    constexpr double h     = 0.001;
    constexpr size_t Niter = 100;

    // 1. determine box from projected/rotated corners

    const std::vector<PointXY> corners{
        {rect.xmin(), rect.ymax()}, {rect.xmax(), rect.ymax()}, {rect.xmax(), rect.ymin()}, {rect.xmin(), rect.ymin()}};

    bool first = true;
    for ( auto& p : corners ) {
        PointLonLat r( p );
        xy2lonlat( r.data() );
        auto rmin = PointLonLat::sub( r, PointLonLat{h, h} );
        auto rmax = PointLonLat::add( r, PointLonLat{h, h} );
        min       = first ? rmin : PointLonLat::componentsMin( min, rmin );
        max       = first ? rmax : PointLonLat::componentsMax( max, rmax );
        first     = false;
    }

    struct crosses_t {
        crosses_t( const ProjectionImpl& projection, const PointXY& p1, const PointXY& p2, const PointXY& p3,
                   const PointXY& p4 ) :
            projection_( projection ),
            quadrilateral_( xy_to_xyz( projection, p1 ), xy_to_xyz( projection, p2 ), xy_to_xyz( projection, p3 ),
                            xy_to_xyz( projection, p4 ) ) {}

        bool operator()( const PointXY& p ) {
            interpolation::method::Ray r( xy_to_xyz( projection_, p ) );
            return quadrilateral_.intersects( r );
        }

    private:
        const ProjectionImpl& projection_;
        interpolation::element::Quad3D quadrilateral_;

        static PointXYZ xy_to_xyz( const ProjectionImpl& projection, const PointXY& p ) {
            PointLonLat q( p );
            projection.xy2lonlat( q.data() );

            PointXYZ r;
            util::Earth::convertSphericalToCartesian( q, r );
            return r;
        }
    } crosses( *this, corners[0], corners[1], corners[2], corners[3] );


    // 2. locate latitude extrema by checking if poles are included (in the unrotated/un-projected frame) and if not,
    // find extrema not at the corners by refining iteratively

    PointXY NP{0, 90};
    PointXY SP{0, -90};
    lonlat2xy( NP.data() );  // unrotate/un-project
    lonlat2xy( SP.data() );

    bool includesNorthPole = crosses( NP );
    bool includesSouthPole = crosses( SP );

    for ( size_t i = 0; i < corners.size(); ++i ) {
        if ( !includesNorthPole || !includesSouthPole ) {
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
                PointLonLat middle( PointXY::middle( A, B ) );
                xy2lonlat( middle.data() );

                min = PointLonLat::componentsMin( min, PointLonLat::sub( middle, PointLonLat{0, h} ) );
                if ( includesSouthPole || is_approximately_lesser_or_equal( min[1], -90. ) ) {
                    includesSouthPole = true;
                    min[1]            = -90.;
                }

                max = PointLonLat::componentsMax( max, PointLonLat::add( middle, PointLonLat{0, h} ) );
                if ( includesNorthPole || is_approximately_lesser_or_equal( 90., max[1] ) ) {
                    includesNorthPole = true;
                    max[1]            = 90.;
                }
            }
        }
    }
    ATLAS_ASSERT( min < max );


    // 3. locate latitude extrema by checking if date line is crossed (in the
    // projected/rotated frame), in which case we assume periodicity and if not, find
    // extrema not at the corners by refining iteratively

    bool crossesDateLine = includesNorthPole || includesSouthPole;

    for ( size_t i = 0; i < corners.size(); ++i ) {
        if ( !crossesDateLine ) {
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
                PointLonLat middle( PointXY::middle( A, B ) );
                xy2lonlat( middle.data() );

                min = PointLonLat::componentsMin( min, PointLonLat::sub( middle, PointLonLat{h, 0} ) );
                max = PointLonLat::componentsMax( max, PointLonLat::add( middle, PointLonLat{h, 0} ) );
                if ( is_approximately_lesser_or_equal( 360., max[0] - min[0] ) ) {
                    crossesDateLine = true;
                }
            }
        }
    }

    if ( crossesDateLine ) {
        max[0] = min[0] + 360.;
    }
    ATLAS_ASSERT( min < max );


    // 4. return bounding box
    return crossesDateLine ? ZonalBandDomain( {min[1], max[1]} )
                           : RectangularDomain( {min[0], max[0]}, {min[1], max[1]} );
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
