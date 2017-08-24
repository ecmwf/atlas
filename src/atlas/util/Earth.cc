/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */


#include "atlas/util/Earth.h"

#include <cmath>
#include "eckit/types/FloatCompare.h"
#include "eckit/exception/Exceptions.h"
#include "atlas/util/Constants.h"
#include "atlas/util/Point.h"


namespace atlas {
namespace util {


static inline double between_m180_and_p180(double a) {
    while (a > 180) {
        a -= 360;
    }
    while (a < -180 ) {
        a += 360;
    }
    return a;
}


double Earth::centralAngle(const PointLonLat& p1, const PointLonLat& p2) {

    /*
     * Δσ = atan( ((cos(ϕ2) * sin(Δλ))^2 + (cos(ϕ1) * sin(ϕ2) - sin(ϕ1) * cos(ϕ2) * cos(Δλ))^2) /
     *            (sin(ϕ1) * sin(ϕ2) + cos(ϕ1) * cos(ϕ2) * cos(Δλ)) )
     *
     * @article{doi:10.1179/sre.1975.23.176.88,
     * author = {T. Vincenty},
     * title = {Direct and Inverse Solutions of Geodesics on the Ellipsoid With Application of Nested Equations},
     * journal = {Survey Review},
     * volume = {23},
     * number = {176},
     * pages = {88-93},
     * year = {1975},
     * doi = {10.1179/sre.1975.23.176.88}
     * }
     */

    ASSERT(-90. <= p1.lat() && p1.lat() <= 90.);
    ASSERT(-90. <= p2.lat() && p2.lat() <= 90.);
    double lambda_deg = between_m180_and_p180(p2.lon() - p1.lon());

    const double
            phi1 = p1.lat() * Constants::degreesToRadians(),
            phi2 = p2.lat() * Constants::degreesToRadians(),
            lambda = lambda_deg * Constants::degreesToRadians(),

            cos_phi1 = std::cos(phi1),
            cos_phi2 = std::cos(phi2),
            sin_phi1 = std::sin(phi1),
            sin_phi2 = std::sin(phi2),
            cos_lambda = std::cos(lambda),
            sin_lambda = std::sin(lambda),

            angle = std::atan2(
                std::sqrt(std::pow(cos_phi2 * sin_lambda, 2) + std::pow(cos_phi1 * sin_phi2 - sin_phi1 * cos_phi2 * cos_lambda, 2)),
                sin_phi1 * sin_phi2 + cos_phi1 * cos_phi2 * cos_lambda );

    if (eckit::types::is_approximately_equal(angle, 0.)) {
        return 0.;
    }

    ASSERT(angle > 0.);
    return angle;
}


double Earth::centralAngle(const PointXYZ& p1, const PointXYZ& p2) {

    // Δσ = 2 * asin( chord / 2 )

    const double d2 = PointXYZ::distance2(p1, p2);
    if (eckit::types::is_approximately_equal(d2, 0.)) {
        return 0.;
    }

    const double
            chord = std::sqrt(d2) / radiusInMeters(),
            angle = std::asin(chord * 0.5) * 2.;

    return angle;
}


double Earth::distanceInMeters(const PointLonLat& p1, const PointLonLat& p2) {
    return radiusInMeters() * centralAngle(p1, p2);
}


double Earth::distanceInMeters(const PointXYZ& p1, const PointXYZ& p2) {
    return radiusInMeters() * centralAngle(p1, p2);
}


void Earth::convertGeodeticToGeocentric(const PointLonLat& p1, PointXYZ& p2, const double& height, const double& radius) {
    convertGeodeticToGeocentric(p1, p2, height, radius, radius);
}


void atlas::util::Earth::convertGeodeticToGeocentric(const atlas::PointLonLat& p1, atlas::PointXYZ& p2, const double& height, const double& a, const double& b) {
    ASSERT(a > 0.);
    ASSERT(b > 0.);
    ASSERT(-90. <= p1.lat() && p1.lat() <= 90.);
    double lambda_deg = between_m180_and_p180(p1.lon());

    // See https://en.wikipedia.org/wiki/Reference_ellipsoid#Coordinates
    // better numerical conditioning for both ϕ (poles) and λ (Greenwich/Date Line)
    const double
            phi = p1.lat() * Constants::degreesToRadians(),
            lambda = lambda_deg * Constants::degreesToRadians(),

            sin_phi = std::sin(phi),
            cos_phi = std::sqrt(1. - sin_phi * sin_phi),
            sin_lambda = std::abs(lambda_deg) < 180. ? std::sin(lambda) : 0.,
            cos_lambda = std::abs(lambda_deg) >  90. ? std::cos(lambda) : std::sqrt(1. - sin_lambda * sin_lambda);

    if (a == b) {
        // no eccentricity
        p2.x() = (a + height) * cos_phi * cos_lambda;
        p2.y() = (a + height) * cos_phi * sin_lambda;
        p2.z() = (a + height) * sin_phi;
    } else {
        const double N_phi = a * a / std::sqrt(a * a * cos_phi * cos_phi + b * b * sin_phi * sin_phi);
        p2.x() = (N_phi + height) * cos_phi * cos_lambda;
        p2.y() = (N_phi + height) * cos_phi * sin_lambda;
        p2.z() = (N_phi * (b * b) / (a * a) + height) * sin_phi;
    }
}


}  // namespace util
}  // namespace atlas
