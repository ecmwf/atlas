/*
 * (C) Copyright 1996- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */


#include "LambertConformalConicProjection.h"

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

namespace {

inline double normalise(double a, double minimum, double globe) {
    while (a >= minimum + globe) {
        a -= globe;
    }
    while (a < minimum) {
        a += globe;
    }
    return a;
}

inline double tan_d(double lat_d) {
    return std::tan(util::Constants::degreesToRadians() * (45. + lat_d * 0.5));
}

inline double sin_d(double theta_d) {
    return std::sin(util::Constants::degreesToRadians() * theta_d);
}

inline double cos_d(double theta_d) {
    return std::cos(util::Constants::degreesToRadians() * theta_d);
}

static ProjectionBuilder<LambertConformalConicProjection> register_1(LambertConformalConicProjection::static_type());

}  // namespace

LambertConformalConicProjection::LambertConformalConicProjection(const eckit::Parametrisation& params) {
    ATLAS_ASSERT(params.get("longitude0", lon0_ = 0));
    ATLAS_ASSERT(params.get("latitude0", lat0_ = 0));
    if (not params.get("latitude1", lat1_)) {
        lat1_ = lat0_;
    }
    if (not params.get("latitude2", lat2_)) {
        lat2_ = lat0_;
    }
    params.get("radius", radius_ = util::Earth::radius());

    n_     = eckit::types::is_approximately_equal(lat1_, lat2_)
                 ? sin_d(lat2_)
                 : std::log(cos_d(lat1_) / cos_d(lat2_)) / std::log(tan_d(lat2_) / tan_d(lat1_));
    inv_n_ = 1. / n_;
    sign_  = n_ < 0. ? -1. : 1.;  // n < 0: adjustment for southern hemisphere

    F_    = (cos_d(lat1_) * std::pow(tan_d(lat1_), n_)) / n_;
    rho0_ = sign_ * radius_ * F_ * std::pow(tan_d(lat0_), -n_);
}

void LambertConformalConicProjection::lonlat2xy(double crd[]) const {
    double rho   = radius_ * F_ * std::pow(tan_d(crd[1]), -n_);
    double theta = n_ * normalise(crd[0] - lon0_, -180, 360);

    crd[XX] = rho * sin_d(theta);
    crd[YY] = rho0_ - rho * cos_d(theta);
}

void LambertConformalConicProjection::xy2lonlat(double crd[]) const {
    double x = sign_ * crd[XX];
    double y = rho0_ - sign_ * crd[YY];

    double rho   = sign_ * std::sqrt(x * x + y * y);
    double theta = std::atan2(x, y) * inv_n_;

    crd[LON] = util::Constants::radiansToDegrees() * theta + lon0_;
    crd[LAT] = eckit::types::is_approximately_equal(rho, 0.)
                   ? 90 * sign_
                   : util::Constants::radiansToDegrees() * 2. * std::atan(std::pow(radius_ * F_ / rho, inv_n_)) - 90.;
}

ProjectionImpl::Jacobian LambertConformalConicProjection::jacobian(const PointLonLat& lonlat) const {
    ProjectionImpl::Jacobian jac;

    const double deg2rad = util::Constants::degreesToRadians();

    double lat = lonlat.lat(), lon = lonlat.lon();
    double tanlat = tan_d(lat), coslat = cos_d(45. + lat * 0.5);

    double rho   = radius_ * F_ * std::pow(tanlat, -n_);
    double theta = n_ * normalise(lon - lon0_, -180, 360);

    double costheta = cos_d(theta), sintheta = sin_d(theta);

    // double x = rho * sintheta;
    // double y = rho0_ - rho * costheta;

    auto cpj = [&](const double dlon, const double dlat, double& dx, double& dy) {
        double drho   = deg2rad * radius_ * F_ * (-n_) * std::pow(tanlat, -n_ - 1) / (coslat * coslat) * dlat * 0.5;
        double dtheta = +n_ * dlon;
        dx            = +drho * sintheta + rho * costheta * deg2rad * dtheta;
        dy            = -drho * costheta + rho * sintheta * deg2rad * dtheta;
    };

    const double dlat = 1.0f;
    const double dlon = 1.0f;

    cpj(dlon, 0, jac[0][0], jac[1][0]);
    cpj(0, dlat, jac[0][1], jac[1][1]);

    return jac;
}

LambertConformalConicProjection::Spec LambertConformalConicProjection::spec() const {
    Spec spec;
    spec.set("type", static_type());
    spec.set("longitude0", lon0_);
    spec.set("latitude0", lat0_);
    spec.set("latitude1", lat1_);
    spec.set("latitude2", lat2_);
    if (!eckit::types::is_approximately_equal(radius_, util::Earth::radius())) {
        spec.set("radius", radius_);
    }
    return spec;
}

void LambertConformalConicProjection::hash(eckit::Hash& h) const {
    h.add(static_type());
    h.add(lat1_);
    h.add(lat2_);
    h.add(lat0_);
    h.add(lon0_);
    h.add(radius_);
}

}  // namespace detail
}  // namespace projection
}  // namespace atlas
