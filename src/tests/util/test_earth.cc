/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#define BOOST_TEST_MODULE TestEarth
#include "ecbuild/boost_test_framework.h"
#include "tests/AtlasFixture.h"

#include <cmath>
#include "atlas/util/Earth.h"
#include "atlas/util/Point.h"

using atlas::util::Earth;

namespace atlas {
namespace test {

BOOST_GLOBAL_FIXTURE( AtlasFixture );
#define CHECK_APPROX_EQUAL(x, y) BOOST_CHECK_LT(std::abs((x) - (y)), 1e-9)
#define CHECK_APPROX_EQUAL_EPS(x, y, eps) BOOST_CHECK_LT(std::abs((x) - (y)), (eps))

BOOST_AUTO_TEST_SUITE( test_earth )

const double R = Earth::radiusInMeters();

// -----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE( test_earth_poles )

BOOST_AUTO_TEST_CASE( test_earth_north_pole )
{
  const PointLonLat p1(0., 90.);
  const PointXYZ p2 = Earth::convertGeodeticToGeocentric(p1);

  BOOST_CHECK_EQUAL(p2.x(), 0);
  BOOST_CHECK_EQUAL(p2.y(), 0);
  BOOST_CHECK_EQUAL(p2.z(), R);
}

BOOST_AUTO_TEST_CASE( test_earth_south_pole )
{
  const PointLonLat p1(0., -90.);
  const PointXYZ p2 = Earth::convertGeodeticToGeocentric(p1);

  BOOST_CHECK_EQUAL(p2.x(), 0);
  BOOST_CHECK_EQUAL(p2.y(), 0);
  BOOST_CHECK_EQUAL(p2.z(), -R);
}

BOOST_AUTO_TEST_SUITE_END()

// -----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE( test_earth_quadrants )

BOOST_AUTO_TEST_CASE( test_earth_lon_0 )
{
  const PointLonLat p1[2] = {{0., 0.}, {-360., 0.}};
  const PointXYZ p2[] = {
      Earth::convertGeodeticToGeocentric(p1[0]),
      Earth::convertGeodeticToGeocentric(p1[1])
  };

  BOOST_CHECK_EQUAL(p2[0].x(), R);
  BOOST_CHECK_EQUAL(p2[0].y(), 0);
  BOOST_CHECK_EQUAL(p2[0].z(), 0);

  BOOST_CHECK(PointXYZ::equal(p2[0], p2[1]));
}

BOOST_AUTO_TEST_CASE( test_earth_lon_90 )
{
  const PointLonLat p1[2] = {{90., 0.}, {-270., 0.}};
  const PointXYZ p2[] = {
      Earth::convertGeodeticToGeocentric(p1[0]),
      Earth::convertGeodeticToGeocentric(p1[1])
  };

  BOOST_CHECK_EQUAL(p2[0].x(), 0);
  BOOST_CHECK_EQUAL(p2[0].y(), R);
  BOOST_CHECK_EQUAL(p2[0].z(), 0);

  BOOST_CHECK(PointXYZ::equal(p2[0], p2[1]));
}

BOOST_AUTO_TEST_CASE( test_earth_lon_180 )
{
  const PointLonLat p1[2] = {{180., 0.}, {-180., 0.}};
  const PointXYZ p2[] = {
      Earth::convertGeodeticToGeocentric(p1[0]),
      Earth::convertGeodeticToGeocentric(p1[1])
  };

  BOOST_CHECK_EQUAL(p2[0].x(), -R);
  BOOST_CHECK_EQUAL(p2[0].y(), 0);
  BOOST_CHECK_EQUAL(p2[0].z(), 0);

  BOOST_CHECK(PointXYZ::equal(p2[0], p2[1]));
}

BOOST_AUTO_TEST_CASE( test_earth_lon_270 )
{
  const PointLonLat p1[2] = {{270., 0.}, {-90., 0.}};
  const PointXYZ p2[] = {
      Earth::convertGeodeticToGeocentric(p1[0]),
      Earth::convertGeodeticToGeocentric(p1[1])
  };

  BOOST_CHECK_EQUAL(p2[0].x(), 0);
  BOOST_CHECK_EQUAL(p2[0].y(), -R);
  BOOST_CHECK_EQUAL(p2[0].z(), 0);

  BOOST_CHECK(PointXYZ::equal(p2[0], p2[1]));
}

BOOST_AUTO_TEST_SUITE_END()

// -----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE( test_earth_octants )

const double L = R * std::sqrt(2) / 2.;

BOOST_AUTO_TEST_CASE( test_earth_lon_45 )
{
  const PointLonLat p1[2] = {{45., 0.}, {-315., 0.}};
  const PointXYZ p2[] = {
      Earth::convertGeodeticToGeocentric(p1[0]),
      Earth::convertGeodeticToGeocentric(p1[1])
  };

  CHECK_APPROX_EQUAL(p2[0].x(), L);
  CHECK_APPROX_EQUAL(p2[0].y(), L);
  BOOST_CHECK_EQUAL(p2[0].z(), 0);

  BOOST_CHECK(PointXYZ::equal(p2[0], p2[1]));
}

BOOST_AUTO_TEST_CASE( test_earth_lon_135 )
{
  const PointLonLat p1[2] = {{135., 0.}, {-225., 0.}};
  const PointXYZ p2[] = {
      Earth::convertGeodeticToGeocentric(p1[0]),
      Earth::convertGeodeticToGeocentric(p1[1])
  };

  CHECK_APPROX_EQUAL(p2[0].x(), -L);
  CHECK_APPROX_EQUAL(p2[0].y(), L);
  BOOST_CHECK_EQUAL(p2[0].z(), 0);

  BOOST_CHECK(PointXYZ::equal(p2[0], p2[1]));
}

BOOST_AUTO_TEST_CASE( test_earth_lon_225 )
{
  const PointLonLat p1[2] = {{225., 0.}, {-135., 0.}};
  const PointXYZ p2[] = {
      Earth::convertGeodeticToGeocentric(p1[0]),
      Earth::convertGeodeticToGeocentric(p1[1])
  };

  CHECK_APPROX_EQUAL(p2[0].x(), -L);
  CHECK_APPROX_EQUAL(p2[0].y(), -L);
  BOOST_CHECK_EQUAL(p2[0].z(), 0);

  BOOST_CHECK(PointXYZ::equal(p2[0], p2[1]));
}

BOOST_AUTO_TEST_CASE( test_earth_lon_315 )
{
  const PointLonLat p1[2] = {{315., 0.}, {-45., 0.}};
  const PointXYZ p2[] = {
      Earth::convertGeodeticToGeocentric(p1[0]),
      Earth::convertGeodeticToGeocentric(p1[1])
  };

  CHECK_APPROX_EQUAL(p2[0].x(), L);
  CHECK_APPROX_EQUAL(p2[0].y(), -L);
  BOOST_CHECK_EQUAL(p2[0].z(), 0);

  BOOST_CHECK(PointXYZ::equal(p2[0], p2[1]));
}

BOOST_AUTO_TEST_SUITE_END()

// -----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE( test_earth_great_circle_navigation )

BOOST_AUTO_TEST_CASE( test_earth_azimuth )
{
    const PointLonLat
            u( 0.,  0.),
            v(90.,  0.),
            w( 0., 90.);

    const double
            a = Earth::azimuth(w, u, v),
            b = Earth::azimuth(v, w, u),
            c = Earth::azimuth(u, v, w);

    CHECK_APPROX_EQUAL(a, M_PI_2);
    CHECK_APPROX_EQUAL(b, M_PI_2);
    CHECK_APPROX_EQUAL(c, M_PI_2);
}

BOOST_AUTO_TEST_CASE( test_earth_great_circle_navigation_article_example )
{
    using namespace std;
    using util::Constants;

    const PointLonLat P1(-71.6, -33.);   // Valparaíso
    const PointLonLat P2(121.8,  31.4);  // Shangai

    const double
            eps = 0.01,
            phi1 = Constants::degreesToRadians() * P1.lat(),
            phi2 = Constants::degreesToRadians() * P2.lat(),

            alpha1  = Earth::course(P1, P2),
            alpha2  = Earth::course(P2, P1) - M_PI,
            sigma12 = Earth::centralAngle(P1, P2),

            sigma01 = atan2(tan(phi1), cos(alpha1)),
            sigma02 = atan2(tan(phi2), cos(alpha2));



    CHECK_APPROX_EQUAL_EPS(Constants::radiansToDegrees() * alpha1,  -94.41, eps);
    CHECK_APPROX_EQUAL_EPS(Constants::radiansToDegrees() * alpha2,  -78.42, eps);
    CHECK_APPROX_EQUAL_EPS(Constants::radiansToDegrees() * sigma12, 168.56, eps);

    CHECK_APPROX_EQUAL_EPS(Constants::radiansToDegrees() * sigma01, -96.76, eps);
    CHECK_APPROX_EQUAL_EPS(Constants::radiansToDegrees() * sigma02,  71.8,  eps);

    double lambda0, alpha0;
    const double
            tan_alpha0 = sin(alpha1) * cos(phi1)
            / sqrt(cos(alpha1) * cos(alpha1) + sin(alpha1) * sin(alpha1) * sin(phi1) * sin(phi1)),
            alpha0____  = atan(tan_alpha0);

    Earth::greatCircleIntersectionAtEquator(P1, P2, lambda0, alpha0);


    CHECK_APPROX_EQUAL_EPS(Constants::radiansToDegrees() * alpha0,   -56.74, eps);
    CHECK_APPROX_EQUAL_EPS(Constants::radiansToDegrees() * lambda0, -169.67, eps);



    // mid-point
    const double lon_rad = Constants::degreesToRadians() * (P1.lon() + P2.lon()) / 2.;
    const double lat_deg = util::Earth::greatCircleLatitudeGivenLongitude(P1, P2, lon_rad);

//    CHECK_APPROX_EQUAL_EPS(lat_deg, 0., eps);
}

BOOST_AUTO_TEST_SUITE_END()

// -----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END()
#undef CHECK_APPROX_EQUAL_EPS
#undef CHECK_APPROX_EQUAL

} // namespace test
} // namespace atlas
