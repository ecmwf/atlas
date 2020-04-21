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
    ProjectionImpl(), rotation_( config ) {}

template <>
void LonLatProjectionT<NotRotated>::xy2lonlat( double[] ) const {}

template <>
void LonLatProjectionT<NotRotated>::lonlat2xy( double[] ) const {}

template <>
RectangularLonLatDomain LonLatProjectionT<NotRotated>::lonlatBoundingBox( const Domain& domain ) const {
    return domain;
}

template <typename Rotation>
RectangularLonLatDomain LonLatProjectionT<Rotation>::lonlatBoundingBox( const Domain& domain ) const {
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

    PointXY NP{xy( {0., 90.} )};
    PointXY SP{xy( {0., -90.} )};

    bounds.includesNorthPole( rect.contains( NP ) );
    bounds.includesSouthPole( rect.contains( SP ) );

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


    // 3. locate latitude extrema by checking if date line is crossed (in the un-projected frame), in which case we
    // assume periodicity and if not, find extrema not at the corners by refining iteratively

    if ( !bounds.crossesDateLine() ) {
        PointLonLat A{xy( {180., -10.} )};
        PointLonLat B{xy( {180., 10.} )};
        eckit::geometry::GreatCircle DL( A, B );

        for ( auto lon : {rect.xmin(), rect.xmax()} ) {
            if ( !bounds.crossesDateLine() ) {
                for ( auto lat : DL.latitude( lon ) ) {
                    if ( ( bounds.crossesDateLine( domain.contains( lon, lat ) ) ) ) {
                        break;
                    }
                }
            }
        }

        for ( auto lat : {rect.ymin(), rect.ymax()} ) {
            if ( !bounds.crossesDateLine() ) {
                for ( auto lon : DL.longitude( lat ) ) {
                    if ( ( bounds.crossesDateLine( domain.contains( lon, lat ) ) ) ) {
                        break;
                    }
                }
            }
        }
    }

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
}  // namespace

// namespace

}  // namespace detail
}  // namespace projection
}  // namespace atlas
