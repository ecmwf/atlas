/*
 * (C) Copyright 1996- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */


#include "LambertAzimuthalEqualAreaProjection.h"

#include <cmath>

#include "eckit/config/Parametrisation.h"
#include "eckit/types/FloatCompare.h"
#include "eckit/utils/Hash.h"

#include "atlas/projection/detail/ProjectionFactory.h"
#include "atlas/runtime/Exception.h"
#include "atlas/util/Config.h"
#include "atlas/util/Constants.h"
#include "atlas/util/CoordinateEnums.h"
#include "atlas/util/Earth.h"

namespace atlas {
namespace projection {
namespace detail {


LambertAzimuthalEqualAreaProjection::LambertAzimuthalEqualAreaProjection( const eckit::Parametrisation& params ) :
    radius_( util::Earth::radius() ) {
    ATLAS_ASSERT( params.get( "central_longitude", reference_[LON] ) );
    ATLAS_ASSERT( params.get( "standard_parallel", reference_[LAT] ) );
    params.get( "radius", radius_ = util::Earth::radius() );
    params.get( "false_northing", false_northing_ );
    params.get( "false_easting", false_easting_ );


    lambda0_  = util::Constants::degreesToRadians() * reference_[LON];
    phi1_     = util::Constants::degreesToRadians() * reference_[LAT];
    sin_phi1_ = std::sin( phi1_ );
    cos_phi1_ = std::cos( phi1_ );
}


void LambertAzimuthalEqualAreaProjection::lonlat2xy( double crd[] ) const {
    double dlambda     = util::Constants::degreesToRadians() * ( crd[LON] - reference_.lon() );
    double cos_dlambda = std::cos( dlambda );
    double sin_dlambda = std::sin( dlambda );

    double phi     = util::Constants::degreesToRadians() * crd[LAT];
    double sin_phi = std::sin( phi );
    double cos_phi = std::cos( phi );

    double kp = radius_ * std::sqrt( 2. / ( 1. + sin_phi1_ * sin_phi + cos_phi1_ * cos_phi * cos_dlambda ) );

    crd[XX] = kp * cos_phi * sin_dlambda;
    crd[YY] = kp * ( cos_phi1_ * sin_phi - sin_phi1_ * cos_phi * cos_dlambda );
    crd[XX] += false_easting_;
    crd[YY] += false_northing_;
}


void LambertAzimuthalEqualAreaProjection::xy2lonlat( double crd[] ) const {
    const double x = crd[XX] - false_easting_;
    const double y = crd[YY] - false_northing_;

    const double rho = std::sqrt( x * x + y * y );
    if ( std::abs( rho ) < 1.e-12 ) {
        crd[LON] = reference_[LON];
        crd[LAT] = reference_[LAT];
        return;
    }

    double c     = 2. * std::asin( rho / ( 2. * radius_ ) );
    double cos_c = std::cos( c );
    double sin_c = std::sin( c );

    double lon_r     = lambda0_ + std::atan2( x * sin_c, rho * cos_phi1_ * cos_c - y * sin_phi1_ * sin_c );
    double sin_lat_r = cos_c * sin_phi1_ + y * sin_c * cos_phi1_ / rho;

    crd[LON] = lon_r * util::Constants::radiansToDegrees();
    crd[LAT] = sin_lat_r > 1 ? 90 : sin_lat_r < -1 ? -90 : util::Constants::radiansToDegrees() * std::asin( sin_lat_r );
}


LambertAzimuthalEqualAreaProjection::Spec LambertAzimuthalEqualAreaProjection::spec() const {
    Spec proj;
    proj.set( "type", static_type() );
    proj.set( "central_longitude", reference_[LON] );
    proj.set( "standard_parallel", reference_[LAT] );
    proj.set( "radius", radius_ );
    proj.set( "false_easting", false_easting_ );
    proj.set( "false_northing", false_northing_ );

    return proj;
}


void LambertAzimuthalEqualAreaProjection::hash( eckit::Hash& h ) const {
    h.add( static_type() );
    h.add( radius_ );
    h.add( reference_[LON] );
    h.add( reference_[LAT] );
}


namespace {
static ProjectionBuilder<LambertAzimuthalEqualAreaProjection> register_1(
    LambertAzimuthalEqualAreaProjection::static_type() );
}


}  // namespace detail
}  // namespace projection
}  // namespace atlas
