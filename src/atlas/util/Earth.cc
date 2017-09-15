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


double Sphere::course(const PointLonLat& p1, const PointLonLat& p2) {
    using namespace std;

    // Solve spherical version of the inverse geodesic problem, see:
    //  https://en.wikipedia.org/wiki/Great-circle_navigation#Course_and_distance
    const double
            lambda12 = Constants::degreesToRadians() * between_m180_and_p180(p2.lon() - p1.lon()),
            phi1     = Constants::degreesToRadians() * p1.lat(),
            phi2     = Constants::degreesToRadians() * p2.lat(),
            a1 = atan2( sin(lambda12),
                        cos(phi1) * tan(phi2) - sin(phi1) * cos(lambda12) );

    return a1;
}


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
            phi1   = Constants::degreesToRadians() * p1.lat(),
            phi2   = Constants::degreesToRadians() * p2.lat(),
            lambda = Constants::degreesToRadians() * lambda_deg,

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


void Sphere::greatCircleIntersectionAtEquator(const PointLonLat& p1, const PointLonLat& p2, double& lambda0, double& alpha0) {
    using namespace std;

    // Estimate great cicle intersection at the equator (lambda0) then solve direct geodesic problem, see:
    //   https://en.wikipedia.org/wiki/Great-circle_navigation#Finding_way-points
    const double
            lambda1 = Constants::degreesToRadians() * p1.lon(),
            phi1    = Constants::degreesToRadians() * p1.lat(),

            tan_phi1 = tan(phi1),
            alpha1 = Sphere::course(p1, p2);

    const double
            cos_alpha1 = cos(alpha1),
            sin_alpha1 = sin(alpha1),
            tan_alpha0 = (sin_alpha1 * cos(phi1))
                   / sqrt(cos_alpha1 * cos_alpha1 + sin_alpha1 * sin_alpha1 * sin(phi1) * sin(phi1)),
            tan_sigma01 = atan(tan_phi1 / cos_alpha1);

    alpha0  = atan(tan_alpha0);
    lambda0 = lambda1 - atan(sin(alpha0) * tan_sigma01);
}


double Sphere::greatCircleLatitudeGivenLongitude(const PointLonLat& p1, const PointLonLat& p2, const double& longitude) {
    using namespace std;

    // Estimate great cicle intersection at the equator (lambda0) then solve direct geodesic problem
    const double
            lambda  = Constants::degreesToRadians() * longitude,
            lambda1 = Constants::degreesToRadians() * p1.lon(),
            phi1    = Constants::degreesToRadians() * p1.lat();

    const double a1 = course(p1, p2);
    if (eckit::types::is_approximately_equal(a1, M_PI_2)) {
        return atan(sin(lambda - lambda1) * tan(phi1));
    }

    double lambda0;
    double a0;
    greatCircleIntersectionAtEquator(p1, p2, lambda0, a0);

    const double
            tan_sigma = tan(lambda - lambda0) / sin(a0),
            sigma = atan(tan_sigma),
            sin_phi = cos(a0) * sin(sigma);

    const double latitude = Constants::radiansToDegrees() * (M_PI_2 - asin(sin_phi));

    static int i = 1;
    std::cout << "h.append(ax.plot("
                 "[" << p1.lon() << ", " << longitude << ", " << p2.lon() << "], "
                 "[" << p1.lat() << ", " << latitude  << ", " << p2.lat() << "], "
                 "'-o', label='Check " << i << "')[0])"
              << std::endl;
    ++i;


    return latitude;
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


double Earth::course(const PointLonLat& source, const PointLonLat& target) {
    return Sphere::course(source, target);
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


void Earth::greatCircleIntersectionAtEquator(const PointLonLat& p1, const PointLonLat& p2, double& lambda0, double& alpha0) {
    Sphere::greatCircleIntersectionAtEquator(p1, p2, lambda0, alpha0);
}


double Earth::greatCircleLatitudeGivenLongitude(const PointLonLat& p1, const PointLonLat& p2, const double& longitude) {
    return Sphere::greatCircleLatitudeGivenLongitude(p1, p2, longitude);
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
