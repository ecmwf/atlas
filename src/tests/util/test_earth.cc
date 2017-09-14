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
#define CHECK_APPROX_EQUAL(x, y) BOOST_CHECK_LT(std::abs((x) - (y)), 1e-9)

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

#undef CHECK_APPROX_EQUAL

BOOST_AUTO_TEST_SUITE_END()

// -----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END()

} // namespace test
} // namespace atlas
