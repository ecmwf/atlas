/*
 * (C) Copyright 1996- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */


#include "atlas/projection/detail/LambertAzimuthalEqualAreaProjection.h"

#include <cmath>

#include "eckit/config/Parametrisation.h"
#include "eckit/exception/Exceptions.h"
#include "eckit/types/FloatCompare.h"
#include "eckit/utils/Hash.h"

#include "atlas/projection/detail/ProjectionFactory.h"
#include "atlas/util/Config.h"
#include "atlas/util/Constants.h"
#include "atlas/util/Earth.h"


namespace atlas {
namespace projection {
namespace detail {


LambertAzimuthalEqualAreaProjection::LambertAzimuthalEqualAreaProjection( const eckit::Parametrisation& params )
    : radius_(util::Earth::radius()) {

    params.get("radius", radius_);
    ASSERT(params.get("central_longitude", reference_[0]));
    ASSERT(params.get("standard_parallel", reference_[1]));

    lambda0_ = util::Constants::degreesToRadians() * reference_[0];
    phi1_    = util::Constants::degreesToRadians() * reference_[1];
    sin_phi1_ = std::sin(phi1_);
    cos_phi1_ = std::cos(phi1_);
}


void LambertAzimuthalEqualAreaProjection::lonlat2xy( double crd[] ) const {
    NOTIMP;
}


void LambertAzimuthalEqualAreaProjection::xy2lonlat( double crd[] ) const {
    const double& x = crd[0];
    const double& y = crd[1];

    const double rho = std::sqrt(x * x + y * y);
    if (eckit::types::is_approximately_equal(rho, 0.)) {
        crd[0] = reference_[0];
        crd[1] = reference_[1];
        return;
    }

    double c = 2. * std::asin(rho / (2. * radius_));
    double cos_c = std::cos(c);
    double sin_c = std::sin(c);

    double lon_r = lambda0_ + std::atan2(x * sin_c, rho * cos_phi1_ * cos_c - y * sin_phi1_ * sin_c);
    double lat_r = std::asin(cos_c * sin_phi1_ + y * sin_c * cos_phi1_ / rho);

    crd[0] = lon_r * util::Constants::radiansToDegrees();
    crd[1] = lat_r * util::Constants::radiansToDegrees();
}


LambertAzimuthalEqualAreaProjection::Spec LambertAzimuthalEqualAreaProjection::spec() const {
    Spec proj;
    proj.set( "type", static_type() );
    proj.set( "central_longitude", reference_[0] );
    proj.set( "standard_parallel", reference_[1] );
    proj.set( "radius", radius_ );

    return proj;
}


void LambertAzimuthalEqualAreaProjection::hash( eckit::Hash& h ) const {
    h.add( static_type() );
    h.add( radius_ );
    h.add( reference_[0] );
    h.add( reference_[1] );
}


namespace {
static ProjectionBuilder<LambertAzimuthalEqualAreaProjection> register_1( LambertAzimuthalEqualAreaProjection::static_type() );
}


}  // namespace detail
}  // namespace projection
}  // namespace atlas

