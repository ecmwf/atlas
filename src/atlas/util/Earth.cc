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

    double lambda_deg = p2.lon() - p1.lon();
    while (lambda_deg > 180.) {
        lambda_deg -= 360.;
    }
    while (lambda_deg < -180.) {
        lambda_deg += 360.;
    }

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
#if 0
    // tests using eckit::geometry::lonlat_to_3d fail with:
    // error: in "test_earth/test_earth_poles/test_earth_north_pole": check p2.x() == 0 has failed [3.9012526007423962e-10 != 0]
    // error: in "test_earth/test_earth_poles/test_earth_south_pole": check p2.x() == 0 has failed [3.9012526007423962e-10 != 0]
    // error: in "test_earth/test_earth_quadrants/test_earth_lon_0": check PointXYZ::equal(p2[0], p2[1]) has failed
    // error: in "test_earth/test_earth_quadrants/test_earth_lon_90": check p2[0].x() == 0 has failed [3.9012526007423962e-10 != 0]
    // error: in "test_earth/test_earth_quadrants/test_earth_lon_90": check PointXYZ::equal(p2[0], p2[1]) has failed
    // error: in "test_earth/test_earth_quadrants/test_earth_lon_180": check p2[0].y() == 0 has failed [7.8025052014847924e-10 != 0]
    // error: in "test_earth/test_earth_quadrants/test_earth_lon_180": check PointXYZ::equal(p2[0], p2[1]) has failed
    // error: in "test_earth/test_earth_quadrants/test_earth_lon_270": check p2[0].x() == 0 has failed [-1.1703757802227187e-09 != 0]
    // error: in "test_earth/test_earth_quadrants/test_earth_lon_270": check PointXYZ::equal(p2[0], p2[1]) has failed
    // error: in "test_earth/test_earth_octants/test_earth_lon_45": check PointXYZ::equal(p2[0], p2[1]) has failed
    // error: in "test_earth/test_earth_octants/test_earth_lon_135": check PointXYZ::equal(p2[0], p2[1]) has failed
    // error: in "test_earth/test_earth_octants/test_earth_lon_225": check PointXYZ::equal(p2[0], p2[1]) has failed
    // error: in "test_earth/test_earth_octants/test_earth_lon_315": check PointXYZ::equal(p2[0], p2[1]) has failed

    double xyz[3] = {0,};
    eckit::geometry::lonlat_to_3d(p1.lon(), p1.lat(), xyz, radius, height);
    p2.x() = xyz[eckit::geometry::XX];
    p2.y() = xyz[eckit::geometry::YY];
    p2.z() = xyz[eckit::geometry::ZZ];
    return;
#endif
    ASSERT(radius > 0.);
    ASSERT(-90. <= p1.lat() && p1.lat() <= 90.);

    double lambda_deg = p1.lon();
    while (lambda_deg > 180.) {
        lambda_deg -= 360.;
    }
    while (lambda_deg <= -180.) {
        lambda_deg += 360.;
    }

    // See https://en.wikipedia.org/wiki/Reference_ellipsoid#Coordinates
    // better numerical conditioning for both ϕ (poles) and λ (Greenwich/Date Line)
    const bool
            lambda_cos_conditioning = std::abs(lambda_deg) <= 90.,
            lambda_sin_conditioning = eckit::types::is_approximately_equal(lambda_deg, 180.);
    const double
            phi = p1.lat() * Constants::degreesToRadians(),
            lambda = lambda_deg * Constants::degreesToRadians(),

            sin_phi = std::sin(phi),
            cos_phi = std::sqrt(1. - sin_phi * sin_phi),
            sin_lambda = !lambda_sin_conditioning ? std::sin(lambda) : 0.,
            cos_lambda = !lambda_cos_conditioning ? std::cos(lambda) : std::sqrt(1. - sin_lambda * sin_lambda);

#if 0
    // WGS84: first numerical eccentricity squared is e2 = 1 - (b*b)/(a*a) = 6.69437999014e-3
    // FIXME correct a, b
    const double
            a = radius,
            b = radius,
            N_phi = a * a / std::sqrt(a * a * cos_phi * cos_phi + b * b * sin_phi * sin_phi);

    p2.x() = (N_phi + height) * cos_phi * cos_lambda;
    p2.y() = (N_phi + height) * cos_phi * sin_lambda;
    p2.z() = (N_phi * (b * b) / (a * a) + height) * sin_phi;
#else
    // ignore eccentricity by setting b = a
    p2.x() = (radius + height) * cos_phi * cos_lambda;
    p2.y() = (radius + height) * cos_phi * sin_lambda;
    p2.z() = (radius + height) * sin_phi;
#endif
}


}  // namespace util
}  // namespace atlas
