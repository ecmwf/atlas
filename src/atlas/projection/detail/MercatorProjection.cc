/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <cmath>
#include <functional>

#include "eckit/config/Parametrisation.h"
#include "eckit/utils/Hash.h"

#include "atlas/projection/detail/MercatorProjection.h"
#include "atlas/projection/detail/ProjectionFactory.h"
#include "atlas/util/Config.h"
#include "atlas/util/Constants.h"
#include "atlas/util/Earth.h"

/*
Projection formula's for Mercator projection from "Map Projections: A Working
Manual"

The origin of the xy-system is at (lon0,0)

*/

namespace {
static double D2R( const double x ) {
    return atlas::util::Constants::degreesToRadians() * x;
}
static double R2D( const double x ) {
    return atlas::util::Constants::radiansToDegrees() * x;
}
}  // namespace

namespace atlas {
namespace projection {
namespace detail {

// constructors
template <typename Rotation>
MercatorProjectionT<Rotation>::MercatorProjectionT( const eckit::Parametrisation& params ) :
    ProjectionImpl(), rotation_( params ) {
    // check presence of radius
    if ( !params.get( "radius", radius_ ) ) {
        radius_ = util::Earth::radius();
    }
    // check presence of lon0
    if ( !params.get( "longitude0", lon0_ ) ) {
        lon0_ = 0.0;
    }

    inv_radius_ = 1. / radius_;
}

template <typename Rotation>
void MercatorProjectionT<Rotation>::lonlat2xy( double crd[] ) const {
    // first unrotate
    rotation_.unrotate( crd );

    // then project
    crd[0] = radius_ * ( D2R( crd[0] - lon0_ ) );
    crd[1] = radius_ * std::log( std::tan( D2R( 45. + crd[1] * 0.5 ) ) );
}

template <typename Rotation>
void MercatorProjectionT<Rotation>::xy2lonlat( double crd[] ) const {
    // first projection
    crd[0] = lon0_ + R2D( crd[0] * inv_radius_ );
    crd[1] = 2. * R2D( std::atan( std::exp( crd[1] * inv_radius_ ) ) ) - 90.;

    // then rotate
    rotation_.rotate( crd );
}

// specification
template <typename Rotation>
typename MercatorProjectionT<Rotation>::Spec MercatorProjectionT<Rotation>::spec() const {
    Spec proj_spec;
    proj_spec.set( "type", static_type() );
    proj_spec.set( "longitude0", lon0_ );
    if ( std::not_equal_to<double>()( radius_, util::Earth::radius() ) ) {
        proj_spec.set( "radius", radius_ );
    }
    rotation_.spec( proj_spec );
    return proj_spec;
}

template <typename Rotation>
void MercatorProjectionT<Rotation>::hash( eckit::Hash& hsh ) const {
    hsh.add( static_type() );
    rotation_.hash( hsh );
    hsh.add( lon0_ );
    hsh.add( radius_ );
}

template class MercatorProjectionT<NotRotated>;
template class MercatorProjectionT<Rotated>;

namespace {
static ProjectionBuilder<MercatorProjection> register_1( MercatorProjection::static_type() );
static ProjectionBuilder<RotatedMercatorProjection> register_2( RotatedMercatorProjection::static_type() );
}  // namespace

}  // namespace detail
}  // namespace projection
}  // namespace atlas
