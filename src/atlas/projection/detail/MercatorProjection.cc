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
#include <limits>
#include <sstream>

#include "eckit/config/Parametrisation.h"
#include "eckit/types/FloatCompare.h"
#include "eckit/utils/Hash.h"

#include "atlas/projection/detail/MercatorProjection.h"
#include "atlas/projection/detail/ProjectionFactory.h"
#include "atlas/runtime/Exception.h"
#include "atlas/util/Config.h"
#include "atlas/util/Constants.h"
#include "atlas/util/CoordinateEnums.h"
#include "atlas/util/Earth.h"
#include "atlas/util/NormaliseLongitude.h"

/*
Projection formula's for Mercator projection from "Map Projections: A Working
Manual"

The origin of the xy-system is at (lon0,0)

*/

namespace {
static constexpr double D2R(const double x) {
    return atlas::util::Constants::degreesToRadians() * x;
}
static constexpr double R2D(const double x) {
    return atlas::util::Constants::radiansToDegrees() * x;
}
}  // namespace

namespace atlas {
namespace projection {
namespace detail {

// constructors
template <typename Rotation>
MercatorProjectionT<Rotation>::MercatorProjectionT(const eckit::Parametrisation& params):
    ProjectionImpl(), normalise_(params), rotation_(params) {
    bool radius_provided = params.get("radius", radius_ = util::Earth::radius());
    k_radius_            = radius_;

    params.get("longitude0", lon0_ = 0.0);
    normalise_mercator_ = util::NormaliseLongitude(lon0_ - 180., lon0_ + 180.);

    if (params.get("latitude1", lat1_ = 0.0)) {
        k_radius_ *= std::cos(D2R(lat1_));
    }

    params.get("false_northing", false_northing_ = 0.);
    params.get("false_easting", false_easting_ = 0.);

    eccentricity_ = 0.;
    auto squared  = [](double x) { return x * x; };
    if (params.get("semi_major_axis", semi_major_axis_ = radius_) &&
        params.get("semi_minor_axis", semi_minor_axis_ = radius_)) {
        ATLAS_ASSERT(
            not radius_provided,
            "Ambiguous parameters provided to MercatorProjection: {radius} and {semi_major_axis,semi_minor_axis}");
        eccentricity_ = std::sqrt(1. - squared(semi_minor_axis_ / semi_major_axis_));
    }

    if (eccentricity_ != 0.) {
        k_radius_ /= std::sqrt(1. - squared(eccentricity_ * std::sin(D2R(lat1_))));
    }

    inv_k_radius_ = 1. / k_radius_;
}

template <typename Rotation>
void MercatorProjectionT<Rotation>::lonlat2xy(double crd[]) const {
    auto t = [&](double& lat) -> double {
        double sinlat = std::sin(D2R(lat));
        double t      = (1. + sinlat) / (1. - sinlat);
        if (eccentricity_ > 0) {  // --> ellipsoidal correction
            double e        = eccentricity_;
            double e_sinlat = e * sinlat;
            t *= std::pow((1. - e_sinlat) / (1. + e_sinlat), e);
        }
        return t;
    };

    // first unrotate
    rotation_.unrotate(crd);

    // then project
    if (eckit::types::is_approximately_equal<double>(crd[LAT], 90., 1e-3)) {
        crd[XX] = false_easting_;
        crd[YY] = std::numeric_limits<double>::infinity();
        return;
    }

    if (eckit::types::is_approximately_equal<double>(crd[LAT], -90., 1e-3)) {
        crd[XX] = false_easting_;
        crd[YY] = -std::numeric_limits<double>::infinity();
        return;
    }

    crd[XX] = k_radius_ * (D2R(normalise_mercator_(crd[LON]) - lon0_));
    crd[YY] = k_radius_ * 0.5 * std::log(t(crd[LAT]));
    crd[XX] += false_easting_;
    crd[YY] += false_northing_;
}

template <typename Rotation>
void MercatorProjectionT<Rotation>::xy2lonlat(double crd[]) const {
    auto compute_lat = [&](double y) -> double {
        //  deepcode ignore FloatingPointEquals: We want exact comparison
        if (eccentricity_ == 0.) {
            return 90. - 2. * R2D(std::atan(std::exp(-y * inv_k_radius_)));
        }
        else {  // eccentricity > 0 --> ellipsoidal correction
            double e = eccentricity_;
            /* Iterative procedure to compute the latitude for the inverse projection
             * From the book "Map Projections-A Working Manual-John P. Snyder (1987)"
             * Equation (7–9) involves rapidly converging iteration: Calculate t from (15-11)
             * Then, assuming an initial trial phi equal to (pi/2 - 2*arctan t) in the right side of equation (7–9),
             * calculate phi on the left side. Substitute the calculated phi into the right side,
             * calculate a new phi, etc., until phi does not change significantly from the preceding trial value of phi
             */
            constexpr double EPSILON = 1.e-12;
            constexpr int MAX_ITER   = 15;

            const double t      = std::exp(-y * inv_k_radius_);
            const double e_half = 0.5 * e;

            double lat = 90. - 2 * R2D(std::atan(t));
            for (int i = 0; i < MAX_ITER; ++i) {
                double e_sinlat = e * std::sin(D2R(lat));
                double dlat =
                    90. - 2. * R2D(std::atan(t * (std::pow(((1.0 - e_sinlat) / (1.0 + e_sinlat)), e_half)))) - lat;
                lat += dlat;
                if (std::abs(dlat) < EPSILON) {
                    return lat;
                }
            }
            ATLAS_THROW_EXCEPTION("Convergence failed in computing latitude in MercatorProjection");
        }
    };

    const double x = crd[XX] - false_easting_;
    const double y = crd[YY] - false_northing_;

    // first projection
    crd[LON] = lon0_ + R2D(x * inv_k_radius_);
    crd[LAT] = compute_lat(y);

    // then rotate
    rotation_.rotate(crd);

    // then normalise
    normalise_(crd);
}

template <typename Rotation>
ProjectionImpl::Jacobian MercatorProjectionT<Rotation>::jacobian(const PointLonLat&) const {
    throw_NotImplemented("MercatorProjectionT::jacobian", Here());
}

// specification
template <typename Rotation>
typename MercatorProjectionT<Rotation>::Spec MercatorProjectionT<Rotation>::spec() const {
    Spec proj;
    proj.set("type", static_type());
    proj.set("longitude0", lon0_);
    proj.set("latitude1", lat1_);
    proj.set("radius", radius_);
    proj.set("false_easting", false_easting_);
    proj.set("false_northing", false_northing_);
    normalise_.spec(proj);
    rotation_.spec(proj);
    return proj;
}

template <typename Rotation>
void MercatorProjectionT<Rotation>::hash(eckit::Hash& hsh) const {
    hsh.add(static_type());
    rotation_.hash(hsh);
    normalise_.hash(hsh);
    hsh.add(lon0_);
    hsh.add(lat1_);
    hsh.add(radius_);
}

template class MercatorProjectionT<NotRotated>;
template class MercatorProjectionT<Rotated>;

namespace {
static ProjectionBuilder<MercatorProjection> register_1(MercatorProjection::static_type());
static ProjectionBuilder<RotatedMercatorProjection> register_2(RotatedMercatorProjection::static_type());
}  // namespace

}  // namespace detail
}  // namespace projection
}  // namespace atlas
