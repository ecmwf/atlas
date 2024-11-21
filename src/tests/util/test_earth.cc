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
#include <limits>

#include "atlas/util/Earth.h"
#include "atlas/util/Point.h"
#include "atlas/util/Geometry.h"

#include "tests/AtlasTestEnvironment.h"

using atlas::util::Earth;

namespace atlas {
namespace test {

const double R = Earth::radius();

// -----------------------------------------------------------------------------
// test_earth_poles

CASE("test_earth_north_pole") {
    const PointLonLat p1(0., 90.);
    PointXYZ p2;
    Earth::convertSphericalToCartesian(p1, p2);

    EXPECT(p2.x() == 0);
    EXPECT(p2.y() == 0);
    EXPECT(p2.z() == R);
}

CASE("test_earth_south_pole") {
    const PointLonLat p1(0., -90.);
    PointXYZ p2;
    Earth::convertSphericalToCartesian(p1, p2);

    EXPECT(p2.x() == 0);
    EXPECT(p2.y() == 0);
    EXPECT(p2.z() == -R);
}

// -----------------------------------------------------------------------------
// test_earth_quadrants

CASE("test_earth_lon_0") {
    const PointLonLat p1[2] = {{0., 0.}, {-360., 0.}};
    PointXYZ p2[2];
    Earth::convertSphericalToCartesian(p1[0], p2[0]);
    Earth::convertSphericalToCartesian(p1[1], p2[1]);

    EXPECT(p2[0].x() == R);
    EXPECT(p2[0].y() == 0);
    EXPECT(p2[0].z() == 0);

    EXPECT(PointXYZ::equal(p2[0], p2[1]));
}

CASE("test_earth_lon_90") {
    const PointLonLat p1[2] = {{90., 0.}, {-270., 0.}};
    PointXYZ p2[2];
    Earth::convertSphericalToCartesian(p1[0], p2[0]);
    Earth::convertSphericalToCartesian(p1[1], p2[1]);

    EXPECT(p2[0].x() == 0);
    EXPECT(p2[0].y() == R);
    EXPECT(p2[0].z() == 0);

    EXPECT(PointXYZ::equal(p2[0], p2[1]));
}

CASE("test_earth_lon_180") {
    const PointLonLat p1[2] = {{180., 0.}, {-180., 0.}};
    PointXYZ p2[2];
    Earth::convertSphericalToCartesian(p1[0], p2[0]);
    Earth::convertSphericalToCartesian(p1[1], p2[1]);

    EXPECT(p2[0].x() == -R);
    EXPECT(p2[0].y() == 0);
    EXPECT(p2[0].z() == 0);

    EXPECT(PointXYZ::equal(p2[0], p2[1]));
}

CASE("test_earth_lon_270") {
    const PointLonLat p1[2] = {{270., 0.}, {-90., 0.}};
    PointXYZ p2[2];
    Earth::convertSphericalToCartesian(p1[0], p2[0]);
    Earth::convertSphericalToCartesian(p1[1], p2[1]);

    EXPECT(p2[0].x() == 0);
    EXPECT(p2[0].y() == -R);
    EXPECT(p2[0].z() == 0);

    EXPECT(PointXYZ::equal(p2[0], p2[1]));
}

// -----------------------------------------------------------------------------
// test_earth_octants

const double L = R * std::sqrt(2) / 2.;

CASE("test_earth_lon_45") {
    const PointLonLat p1[2] = {{45., 0.}, {-315., 0.}};
    PointXYZ p2[2];
    Earth::convertSphericalToCartesian(p1[0], p2[0]);
    Earth::convertSphericalToCartesian(p1[1], p2[1]);

    EXPECT(eckit::types::is_approximately_equal(p2[0].x(), L));
    EXPECT(eckit::types::is_approximately_equal(p2[0].y(), L));
    EXPECT(p2[0].z() == 0);

    EXPECT(PointXYZ::equal(p2[0], p2[1]));
}

CASE("test_earth_lon_135") {
    const PointLonLat p1[2] = {{135., 0.}, {-225., 0.}};
    PointXYZ p2[2];
    Earth::convertSphericalToCartesian(p1[0], p2[0]);
    Earth::convertSphericalToCartesian(p1[1], p2[1]);

    EXPECT(eckit::types::is_approximately_equal(p2[0].x(), -L));
    EXPECT(eckit::types::is_approximately_equal(p2[0].y(), L));
    EXPECT(p2[0].z() == 0);

    EXPECT(PointXYZ::equal(p2[0], p2[1]));
}

CASE("test_earth_lon_225") {
    const PointLonLat p1[2] = {{225., 0.}, {-135., 0.}};
    PointXYZ p2[2];
    Earth::convertSphericalToCartesian(p1[0], p2[0]);
    Earth::convertSphericalToCartesian(p1[1], p2[1]);

    EXPECT(eckit::types::is_approximately_equal(p2[0].x(), -L));
    EXPECT(eckit::types::is_approximately_equal(p2[0].y(), -L));
    EXPECT(p2[0].z() == 0);

    EXPECT(PointXYZ::equal(p2[0], p2[1]));
}

CASE("test_earth_lon_315") {
    const PointLonLat p1[2] = {{315., 0.}, {-45., 0.}};
    PointXYZ p2[2];
    Earth::convertSphericalToCartesian(p1[0], p2[0]);
    Earth::convertSphericalToCartesian(p1[1], p2[1]);

    EXPECT(eckit::types::is_approximately_equal(p2[0].x(), L));
    EXPECT(eckit::types::is_approximately_equal(p2[0].y(), -L));
    EXPECT(p2[0].z() == 0);

    EXPECT(PointXYZ::equal(p2[0], p2[1]));
}

CASE("test_earth_great_circle_latitude_given_longitude") {
    // latitude at Valparaíso-Shanghai mid-point
    const PointLonLat P1(-71.6, -33.);
    const PointLonLat P2(121.8, 31.4);

    double lon = -159.18;
    double lat = Earth::greatCircleLatitudeGivenLongitude(P1, P2, lon);

    EXPECT(eckit::types::is_approximately_equal(lat, -6.81, 0.01));
}

CASE("test_across_pole") {
    const PointLonLat P1{-16.875,-105.255};
    atlas::geometry::Earth geo;
    EXPECT_APPROX_EQ(geo.xyz(P1), PointXYZ(-1.60418e+06,486624,-6.14673e+06), 50.);
}


}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
