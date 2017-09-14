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


PointXYZ EllipsoidOfRevolution::convertEllipsoidalToCartesian(const atlas::PointLonLat& p, const double& a, const double& b, const double& height) {
    ASSERT(a > 0.);
    ASSERT(b > 0.);
    ASSERT(-90. <= p.lat() && p.lat() <= 90.);

    double lambda_deg = between_m180_and_p180(p.lon());

    // See https://en.wikipedia.org/wiki/Reference_ellipsoid#Coordinates
    // numerical conditioning for both ϕ (poles) and λ (Greenwich/Date Line)
    const double
            phi = p.lat() * Constants::degreesToRadians(),
            lambda = lambda_deg * Constants::degreesToRadians(),

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


//------------------------------------------------------------------------------------------------------


double Sphere::centralAngle(const PointLonLat& p1, const PointLonLat& p2) {
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


double Sphere::greatCircleLatitudeGivenLongitude(const PointLonLat& p1, const PointLonLat& p2, const double& longitude, const double& radius) {
    ASSERT(radius > 0.);
    using namespace std;
    using namespace eckit::types;

    if (is_approximately_equal(p1.lat(), p2.lat())) {
        return p1.lat();
    }

    const double
            lambda   = Constants::degreesToRadians() * longitude,
            lambda1  = Constants::degreesToRadians() * p1.lon(),
            lambda12 = Constants::degreesToRadians() * between_m180_and_p180(p2.lon() - p1.lon()),
            phi1     = Constants::degreesToRadians() * p1.lat(),
            phi2     = Constants::degreesToRadians() * p2.lat();

    // 1. solve spherical triangle using the North pole (two sides and the included angle given, SAS),
    // 2. estimate great cicle intersection at the equator (lambda0),
    // 3. finally solve direct geodesic problem
    // See:
    //   https://en.wikipedia.org/wiki/Solution_of_triangles#Two_sides_and_the_included_angle_given
    //   https://en.wikipedia.org/wiki/Great-circle_navigation#Finding_way-points
    static const PointLonLat NP(0., 90.);

    const double
            sin_phi1 = sin(phi1),
            cos_phi1 = cos(phi1),
            cos_phi2 = cos(phi2),

            a = distanceInMeters(NP, p1, radius),
            b = distanceInMeters(NP, p2, radius),
            c = sqrt(a * a + b * b - 2 * a * b * cos(lambda12)),

            a1 = acos( min(1., max(-1., (b * b + c * c - a * a) / (2. * b * c) ))),
            cos_a1 = cos(a1),
            sin_a1 = sin(a1);

    // don't handle points at the poles nor great circles that are meridians
    ASSERT(!is_approximately_equal(cos_phi1, 0.));
    ASSERT(!is_approximately_equal(cos_phi2, 0.));
    ASSERT(!is_approximately_equal(sin_a1, 0.));

    const double
            tan_s01 = is_approximately_equal(a1, M_PI_2) ? 0. : tan(phi1) / cos(a1),
            tan_a0 = (sin_a1 * cos_phi1)
                   / sqrt(cos_a1 * cos_a1 + sin_a1 * sin_a1 * sin_phi1 * sin_phi1),
            sin_a0 = sin(atan(tan_a0)),

            lambda01 = atan(sin_a0 * tan_s01),
            lambda0  = lambda1 - lambda01,

            tan_phi = sin(lambda - lambda0) / tan_a0,
            phi = atan(tan_phi);


std::cout << "    \"\\n\" \"ax.plot([" << p1.lon() << "," << lambda0 * Constants::radiansToDegrees() << "," << p2.lon() << "], [" << p1.lat() << ",0.," << p2.lat() << "], '-o')\"" << std::endl;



    return Constants::radiansToDegrees() * phi;
}


PointXYZ Sphere::convertSphericalToCartesian(const PointLonLat& p, const double& radius, const double& height) {
    return EllipsoidOfRevolution::convertEllipsoidalToCartesian(p, radius, radius, height);
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


double Earth::greatCircleLatitudeGivenLongitude(const PointLonLat& p1, const PointLonLat& p2, const double& longitude, const double& radius) {
    return Sphere::greatCircleLatitudeGivenLongitude(p1, p2, longitude, radius);
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
