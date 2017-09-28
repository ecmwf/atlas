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


//------------------------------------------------------------------------------------------------------


static inline double between_m180_and_p180(double a) {
    while (a > 180) {
        a -= 360;
    }
    while (a < -180 ) {
        a += 360;
    }
    return a;
}


//------------------------------------------------------------------------------------------------------


double Sphere::azimuth(const PointLonLat& u, const PointLonLat& v, const PointLonLat& w) {
    using namespace std;

    // Solve spherical triangle using law of cosines (cosine rule for sides), see:
    //   https://en.wikipedia.org/wiki/Spherical_law_of_cosines
    const double
            a = distanceInMeters(u, v, 1.),
            b = distanceInMeters(u, w, 1.),
            c = distanceInMeters(v, w, 1.),
            cos_C = - (cos(a) * cos(b) - cos(c)) /  (sin(a) * sin(b));

    return acos( min(1., max( -1., cos_C )));
}


double Sphere::centralAngle(const PointLonLat& p1, const PointLonLat& p2) {
    using namespace std;

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

    const double
            phi1   = Constants::degreesToRadians() * p1.lat(),
            phi2   = Constants::degreesToRadians() * p2.lat(),
            lambda = Constants::degreesToRadians() * (p2.lon() - p1.lon()),

            cos_phi1 = cos(phi1),
            sin_phi1 = sin(phi1),
            cos_phi2 = cos(phi2),
            sin_phi2 = sin(phi2),
            cos_lambda = cos(lambda),
            sin_lambda = sin(lambda),

            angle = atan2(
                sqrt(pow(cos_phi2 * sin_lambda, 2) + pow(cos_phi1 * sin_phi2 - sin_phi1 * cos_phi2 * cos_lambda, 2)),
                sin_phi1 * sin_phi2 + cos_phi1 * cos_phi2 * cos_lambda );

    if (eckit::types::is_approximately_equal(angle, 0.)) {
        return 0.;
    }

    ASSERT(angle > 0.);
    return angle;
}


double Sphere::centralAngle(const PointXYZ& p1, const PointXYZ& p2, const double& radius) {
    ASSERT(radius > 0.);

    // Δσ = 2 * asin( chord / 2 )

    const double d2 = PointXYZ::distance2(p1, p2);
    if (eckit::types::is_approximately_equal(d2, 0.)) {
        return 0.;
    }

    const double
            chord = std::sqrt(d2) / radius,
            angle = std::asin(chord * 0.5) * 2.;

    return angle;
}


double Sphere::distanceInMeters(const PointLonLat& p1, const PointLonLat& p2, const double& radius) {
    return radius * centralAngle(p1, p2);
}


double Sphere::distanceInMeters(const PointXYZ& p1, const PointXYZ& p2, const double& radius) {
    return radius * centralAngle(p1, p2, radius);
}


void Sphere::greatCircleLatitudeGivenLongitude(const PointLonLat& p1, const PointLonLat& p2, PointLonLat& p) {
    using namespace std;

    //  Intermediate great circle points (not applicable for meridians), see
    //   http://www.edwilliams.org/avform.htm#Int
    ASSERT(!eckit::types::is_approximately_equal(p1.lon(), p2.lon()));
    const double
            phi1     = Constants::degreesToRadians() * p1.lat(),
            phi2     = Constants::degreesToRadians() * p2.lat(),
            lambda1p = Constants::degreesToRadians() * (p.lon() - p1.lon()),
            lambda2p = Constants::degreesToRadians() * (p.lon() - p2.lon()),
            lambda   = Constants::degreesToRadians() * (p2.lon() - p1.lon());

    p.lat() = Constants::radiansToDegrees() * atan(
                (tan(phi2) * sin(lambda1p) - tan(phi1) * sin(lambda2p)) /
                (sin(lambda)) );
}


PointXYZ Sphere::convertSphericalToCartesian(const PointLonLat& p, const double& radius, const double& height) {
    ASSERT(radius > 0.);
    ASSERT(-90. <= p.lat() && p.lat() <= 90.);

    // See https://en.wikipedia.org/wiki/Reference_ellipsoid#Coordinates
    // numerical conditioning for both ϕ (poles) and λ (Greenwich/Date Line)
    const double
            &a = radius,
            &b = radius,

            lambda_deg = between_m180_and_p180(p.lon()),
            lambda = Constants::degreesToRadians() * lambda_deg,
            phi    = Constants::degreesToRadians() * p.lat(),

            sin_phi = std::sin(phi),
            cos_phi = std::sqrt(1. - sin_phi * sin_phi),
            sin_lambda = std::abs(lambda_deg) < 180. ? std::sin(lambda) : 0.,
            cos_lambda = std::abs(lambda_deg) >  90. ? std::cos(lambda) : std::sqrt(1. - sin_lambda * sin_lambda);

    if (a == b) {
        // no eccentricity
        return PointXYZ(
                    (a + height) * cos_phi * cos_lambda,
                    (a + height) * cos_phi * sin_lambda,
                    (a + height) * sin_phi );
    }

    const double N_phi = a * a / std::sqrt(a * a * cos_phi * cos_phi + b * b * sin_phi * sin_phi);
    return PointXYZ(
        (N_phi + height) * cos_phi * cos_lambda,
        (N_phi + height) * cos_phi * sin_lambda,
        (N_phi * (b * b) / (a * a) + height) * sin_phi );
}


PointLonLat Sphere::convertCartesianToSpherical(const PointXYZ& p, const double& radius) {
    ASSERT(radius > 0.);

    // numerical conditioning for both z (poles) and y
    const double
            x = p.x(),
            y = eckit::types::is_approximately_equal(p.y(), 0.) ? 0. : p.y(),
            z = std::min(radius, std::max(-radius, p.z())) / radius;

    return PointLonLat(
                Constants::radiansToDegrees() * std::atan2(y, x),
                Constants::radiansToDegrees() * std::asin(z) );
}


//------------------------------------------------------------------------------------------------------


double Earth::azimuth(const PointLonLat& source, const PointLonLat& target, const PointLonLat& reference) {
    return Sphere::azimuth(source, target, reference);
}


double Earth::centralAngle(const PointLonLat& p1, const PointLonLat& p2) {
    return Sphere::centralAngle(p1, p2);
}


double Earth::centralAngle(const PointXYZ& p1, const PointXYZ& p2, const double& radius) {
    return Sphere::centralAngle(p1, p2, radius);
}


double Earth::distanceInMeters(const PointLonLat& p1, const PointLonLat& p2, const double& radius) {
    return Sphere::distanceInMeters(p1, p2, radius);
}


double Earth::distanceInMeters(const PointXYZ& p1, const PointXYZ& p2, const double& radius) {
    return Sphere::distanceInMeters(p1, p2, radius);
}


void Earth::greatCircleLatitudeGivenLongitude(const PointLonLat& p1, const PointLonLat& p2, PointLonLat& p) {
    Sphere::greatCircleLatitudeGivenLongitude(p1, p2, p);
}


PointXYZ Earth::convertGeodeticToGeocentric(const PointLonLat& p, const double& radius, const double& height) {
    return Sphere::convertSphericalToCartesian(p, radius, height);
}


PointLonLat Earth::convertGeocentricToGeodetic(const PointXYZ& p, const double& radius) {
    return Sphere::convertCartesianToSpherical(p, radius);
}


//------------------------------------------------------------------------------------------------------


}  // namespace util
}  // namespace atlas
