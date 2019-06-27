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

#include "eckit/geometry/GreatCircle.h"
#include "eckit/types/FloatCompare.h"
#include "eckit/utils/Hash.h"

#include "atlas/domain.h"
#include "atlas/projection/detail/LonLatProjection.h"
#include "atlas/projection/detail/ProjectionFactory.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"
#include "atlas/util/Config.h"

namespace atlas {
namespace projection {
namespace detail {

namespace {

struct Derivate {
    Derivate( const util::Rotation& rotation, PointLonLat H ) : rotation_( rotation ), H_( H ) {}
    virtual PointLonLat d( PointLonLat ) const = 0;
    const util::Rotation& rotation_;
    const PointLonLat H_;
    void xy2lonlat( double crd[] ) const { rotation_.rotate( crd ); }
};

struct DerivateForwards final : Derivate {
    using Derivate::Derivate;
    PointLonLat d( PointLonLat P ) const override {
        PointXY A( P );
        PointXY B( PointXY::add( P, Derivate::H_ ) );
        Derivate::xy2lonlat( A.data() );
        Derivate::xy2lonlat( B.data() );
        return PointXY::div( PointXY::sub( B, A ), PointXY::norm( Derivate::H_ ) );
    }
};

struct DerivateBackwards final : Derivate {
    using Derivate::Derivate;
    PointLonLat d( PointLonLat P ) const override {
        PointXY A( PointXY::sub( P, Derivate::H_ ) );
        PointXY B( P );
        Derivate::xy2lonlat( A.data() );
        Derivate::xy2lonlat( B.data() );
        return PointXY::div( PointXY::sub( B, A ), PointXY::norm( Derivate::H_ ) );
    }
};

struct DerivateCentral final : Derivate {
    using Derivate::Derivate;
    PointLonLat d( PointLonLat P ) const override {
        PointXY H2( PointXY::mul( Derivate::H_, 0.5 ) );
        PointXY A( PointXY::sub( P, H2 ) );
        PointXY B( PointXY::add( P, H2 ) );
        Derivate::xy2lonlat( A.data() );
        Derivate::xy2lonlat( B.data() );
        return PointXY::div( PointXY::sub( B, A ), PointXY::norm( Derivate::H_ ) );
    }
};

}  // (anonymous namespace)

template <typename Rotation>
LonLatProjectionT<Rotation>::LonLatProjectionT( const eckit::Parametrisation& config ) :
    ProjectionImpl(),
    rotation_( config ) {}

template <>
PointLonLat LonLatProjectionT<NotRotated>::rotate( const PointXY& p ) const {
    return {p};
}

template <>
PointXY LonLatProjectionT<NotRotated>::unrotate( const PointLonLat& q ) const {
    return {q};
}

template <typename Rotation>
PointLonLat LonLatProjectionT<Rotation>::rotate( const PointXY& p ) const {
    PointLonLat q( p );
    this->xy2lonlat( q.data() );
    return q;
}

template <typename Rotation>
PointXY LonLatProjectionT<Rotation>::unrotate( const PointLonLat& q ) const {
    PointLonLat p( q );
    this->lonlat2xy( p.data() );
    return p;
}

template <>
Domain LonLatProjectionT<NotRotated>::boundingBox( const Domain& domain ) const {
    return domain;
}

template <typename Rotation>
Domain LonLatProjectionT<Rotation>::boundingBox( const Domain& domain ) const {
    using eckit::types::is_approximately_lesser_or_equal;
    using eckit::types::is_strictly_greater;

    RectangularDomain rect( domain );
    ATLAS_ASSERT( rect );

    // 0. setup
    Point2 min, max;

    constexpr double h = 0.001;

    // 1. determine box from rotated corners

    const std::vector<PointLonLat> corners{
        {rect.xmin(), rect.ymax()}, {rect.xmax(), rect.ymax()}, {rect.xmax(), rect.ymin()}, {rect.xmin(), rect.ymin()}};

    bool first = true;
    for ( auto& p : corners ) {
        PointLonLat r( rotate( p ) );
        min   = first ? r : Point2::componentsMin( min, r );
        max   = first ? r : Point2::componentsMax( max, r );
        first = false;
    }


    // 2. locate latitude extrema by checking if poles are included (in the
    // unrotated frame) and if not, find extrema not at the corners by refining
    // iteratively

    PointLonLat NP{unrotate( {0., 90.} )};
    PointLonLat SP{unrotate( {0., -90.} )};

    bool includesNorthPole = domain.contains( NP.lon(), NP.lat() );
    bool includesSouthPole = domain.contains( SP.lon(), SP.lat() );

    if ( !includesNorthPole || !includesSouthPole ) {
        for ( size_t i = 0; i < corners.size(); ++i ) {
            PointLonLat A = corners[i];
            PointLonLat B = corners[( i + 1 ) % corners.size()];

            // finite difference vector derivative (H is the perturbation vector)
            const PointLonLat H{Point2::mul( Point2::normalize( Point2::sub( B, A ) ), h )};
            std::unique_ptr<Derivate> derivate( new DerivateCentral( this->rotation_, H) );

            double derivativeAtA = derivate->d( A ).lat();
            double derivativeAtB = derivate->d( B ).lat();

            if ( !is_strictly_greater( derivativeAtA * derivativeAtB, 0. ) ) {
                for ( size_t cnt = 0; cnt < 100; ++cnt ) {
                    PointLonLat M        = PointLonLat::middle( A, B );
                    double derivativeAtM = derivate->d( M ).lat();
                    if ( is_strictly_greater( derivativeAtA * derivativeAtM, 0. ) ) {
                        A             = M;
                        derivativeAtA = derivativeAtM;
                    }
                    else if ( is_strictly_greater( derivativeAtB * derivativeAtM, 0. ) ) {
                        B             = M;
                        derivativeAtB = derivativeAtM;
                    }
                    else {
                        break;
                    }
                }

                PointLonLat r( rotate( PointLonLat::middle( A, B ) ) );
                min = Point2::componentsMin( min, r );
                max = Point2::componentsMax( max, r );
            }
        }

        // extend by 'a small amount' (arbitrary)
        min = Point2::sub( min, Point2{0, h} );
        max = Point2::add( max, Point2{0, h} );

        includesNorthPole = includesNorthPole || is_approximately_lesser_or_equal( 90., max[1] );
        includesSouthPole = includesSouthPole || is_approximately_lesser_or_equal( min[1], -90. );
    }
    ATLAS_ASSERT( min[1] < max[1] );


    // 3. locate latitude extrema by checking if date line is crossed (in the
    // unrotated frame), in which case we assume periodicity and if not, find
    // extrema not at the corners by refining iteratively

    bool crossesDateLine = includesNorthPole || includesSouthPole;

    if ( !crossesDateLine ) {
        PointLonLat A{unrotate( {180., -10} )};
        PointLonLat B{unrotate( {180., 10} )};
        eckit::geometry::GreatCircle DL( A, B );

        for ( auto lon : {rect.xmin(), rect.xmax()} ) {
            if ( !crossesDateLine ) {
                for ( auto lat : DL.latitude( lon ) ) {
                    if ( ( crossesDateLine = domain.contains( lon, lat ) ) ) {
                        break;
                    }
                }
            }
        }

        for ( auto lat : {rect.ymin(), rect.ymax()} ) {
            if ( !crossesDateLine ) {
                for ( auto lon : DL.longitude( lat ) ) {
                    if ( ( crossesDateLine = domain.contains( lon, lat ) ) ) {
                        break;
                    }
                }
            }
        }
    }

    if ( !crossesDateLine ) {
        for ( size_t i = 0; i < corners.size(); ++i ) {
            PointLonLat A = corners[i];
            PointLonLat B = corners[( i + 1 ) % corners.size()];

            // finite difference vector derivative (H is the perturbation vector)
            const PointLonLat H{Point2::mul( Point2::normalize( Point2::sub( B, A ) ), h )};
            std::unique_ptr<Derivate> derivate( new DerivateCentral( this->rotation_, H) );

            double derivativeAtA = derivate->d( A ).lon();
            double derivativeAtB = derivate->d( B ).lon();

            if ( !is_strictly_greater( derivativeAtA * derivativeAtB, 0. ) ) {
                for ( size_t cnt = 0; cnt < 100; ++cnt ) {
                    PointLonLat M        = PointLonLat::middle( A, B );
                    double derivativeAtM = derivate->d( M ).lon();
                    if ( is_strictly_greater( derivativeAtA * derivativeAtM, 0. ) ) {
                        A             = M;
                        derivativeAtA = derivativeAtM;
                    }
                    else if ( is_strictly_greater( derivativeAtB * derivativeAtM, 0. ) ) {
                        B             = M;
                        derivativeAtB = derivativeAtM;
                    }
                    else {
                        break;
                    }
                }

                PointLonLat r( rotate( PointLonLat::middle( A, B ) ) );
                min = Point2::componentsMin( min, r );
                max = Point2::componentsMax( max, r );
            }
        }

        // extend by 'a small amount' (arbitrary)
        min = Point2::sub( min, Point2{h, 0} );
        max = Point2::add( max, Point2{h, 0} );

        crossesDateLine = is_approximately_lesser_or_equal( 360., max[0] - min[0] );
    }
    ATLAS_ASSERT( min[0] < max[0] );


    // 4. set bounding box
    double n = includesNorthPole ? 90. : max[1];
    double s = includesSouthPole ? -90. : min[1];
    double w = crossesDateLine ? 0 : min[0];
    double e = crossesDateLine ? 360. : max[0];

    return RectangularDomain({w, e}, {s, n});
}

template <typename Rotation>
typename LonLatProjectionT<Rotation>::Spec LonLatProjectionT<Rotation>::spec() const {
    Spec proj_spec;
    proj_spec.set( "type", static_type() );
    rotation_.spec( proj_spec );
    return proj_spec;
}

template <typename Rotation>
void LonLatProjectionT<Rotation>::hash( eckit::Hash& hsh ) const {
    hsh.add( static_type() );
    rotation_.hash( hsh );
}

template class LonLatProjectionT<NotRotated>;
template class LonLatProjectionT<Rotated>;

namespace {
static ProjectionBuilder<LonLatProjection> register_1( LonLatProjection::static_type() );
static ProjectionBuilder<RotatedLonLatProjection> register_2( RotatedLonLatProjection::static_type() );
}

// namespace

}  // namespace detail
}  // namespace projection
}  // namespace atlas
