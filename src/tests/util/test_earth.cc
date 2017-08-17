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

#include "atlas/util/Earth.h"
#include "atlas/util/Point.h"

using atlas::util::Earth;

namespace atlas {
namespace test {

BOOST_GLOBAL_FIXTURE( AtlasFixture );

BOOST_AUTO_TEST_SUITE( test_earth )

// -----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( test_earth_lon_0_lat_0 )
{
  const PointLonLat p1(0., 0.);
  PointXYZ p2;
  Earth::convertGeodeticToGeocentric(p1, p2);

  BOOST_CHECK_EQUAL(p2.x(), Earth::radiusInMeters());
  BOOST_CHECK_EQUAL(p2.y(), 0);
  BOOST_CHECK_EQUAL(p2.z(), 0);
}

BOOST_AUTO_TEST_CASE( test_earth_lon_p180_lat_0 )
{
  const PointLonLat p1(180., 0.);
  PointXYZ p2;
  Earth::convertGeodeticToGeocentric(p1, p2);

  BOOST_CHECK_EQUAL(p2.x(), -Earth::radiusInMeters());
  BOOST_CHECK_EQUAL(p2.y(), 0);
  BOOST_CHECK_EQUAL(p2.z(), 0);
}

BOOST_AUTO_TEST_CASE( test_earth_lon_m180_lat_0 )
{
  const PointLonLat p1(-180., 0.);
  PointXYZ p2;
  Earth::convertGeodeticToGeocentric(p1, p2);

  BOOST_CHECK_EQUAL(p2.x(), -Earth::radiusInMeters());
  BOOST_CHECK_EQUAL(p2.y(), 0);
  BOOST_CHECK_EQUAL(p2.z(), 0);
}

BOOST_AUTO_TEST_CASE( test_earth_lon_p90_lat_0 )
{
  const PointLonLat p1(90., 0.);
  PointXYZ p2;
  Earth::convertGeodeticToGeocentric(p1, p2);

  BOOST_CHECK_EQUAL(p2.x(), 0);
  BOOST_CHECK_EQUAL(p2.y(), Earth::radiusInMeters());
  BOOST_CHECK_EQUAL(p2.z(), 0);
}

BOOST_AUTO_TEST_CASE( test_earth_lon_m90_lat_0 )
{
  const PointLonLat p1(-90., 0.);
  PointXYZ p2;
  Earth::convertGeodeticToGeocentric(p1, p2);

  BOOST_CHECK_EQUAL(p2.x(), 0);
  BOOST_CHECK_EQUAL(p2.y(), -Earth::radiusInMeters());
  BOOST_CHECK_EQUAL(p2.z(), 0);
}

BOOST_AUTO_TEST_CASE( test_earth_lon_0_lat_p90 )
{
  const PointLonLat p1(0., 90.);
  PointXYZ p2;
  Earth::convertGeodeticToGeocentric(p1, p2);

  BOOST_CHECK_EQUAL(p2.x(), 0);
  BOOST_CHECK_EQUAL(p2.y(), 0);
  BOOST_CHECK_EQUAL(p2.z(), Earth::radiusInMeters());
}

BOOST_AUTO_TEST_CASE( test_earth_lon_0_lat_m90 )
{
  const PointLonLat p1(0., -90.);
  PointXYZ p2;
  Earth::convertGeodeticToGeocentric(p1, p2);

  BOOST_CHECK_EQUAL(p2.x(), 0);
  BOOST_CHECK_EQUAL(p2.y(), 0);
  BOOST_CHECK_EQUAL(p2.z(), -Earth::radiusInMeters());
}

// -----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END()

} // namespace test
} // namespace atlas
