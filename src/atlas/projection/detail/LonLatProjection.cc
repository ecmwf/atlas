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

    // 1. determine box from rotated corners

    const std::vector<PointXY> corners{
        {rect.xmin(), rect.ymax()}, {rect.xmax(), rect.ymax()}, {rect.xmax(), rect.ymin()}, {rect.xmin(), rect.ymin()}};

    bool first = true;
    for ( auto& p : corners ) {
        PointLonLat r( rotate( p ) );
        auto rmin = PointLonLat::sub( r, PointLonLat{h, h} );
        auto rmax = PointLonLat::add( r, PointLonLat{h, h} );
        min       = first ? rmin : PointLonLat::componentsMin( min, rmin );
        max       = first ? rmax : PointLonLat::componentsMax( max, rmax );
        first     = false;
    }

    // 2. locate latitude extrema by checking if poles are included (in the unrotated frame) and if not, find extrema
    // not at the corners by refining iteratively

    PointXY NP{unrotate( {0., 90.} )};
    PointXY SP{unrotate( {0., -90.} )};

    bool includesNorthPole = rect.contains( NP );
    bool includesSouthPole = rect.contains( SP );

    for ( size_t i = 0; i < corners.size(); ++i ) {
        if ( !includesNorthPole || !includesSouthPole ) {
            PointXY A = corners[i];
            PointXY B = corners[( i + 1 ) % corners.size()];

            // finite difference vector derivative (H is the perturbation vector)
            const PointXY H{PointXY::mul( PointXY::normalize( PointXY::sub( B, A ) ), h )};
            std::unique_ptr<Derivate> derivate( DerivateFactory::build( "central", *this, H ) );

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
                PointLonLat middle( rotate( PointXY::middle( A, B ) ) );

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

    for ( size_t i = 0; i < corners.size(); ++i ) {
        if ( !crossesDateLine ) {
            PointXY A = corners[i];
            PointXY B = corners[( i + 1 ) % corners.size()];

            // finite difference vector derivative (H is the perturbation vector)
            const PointXY H{PointXY::mul( PointXY::normalize( PointXY::sub( B, A ) ), h )};
            std::unique_ptr<Derivate> derivate( DerivateFactory::build( "central", *this, H ) );

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
                PointLonLat middle( rotate( PointXY::middle( A, B ) ) );

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
