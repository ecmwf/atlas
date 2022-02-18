/*
 * (C) Crown Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "eckit/types/FloatCompare.h"

#include "atlas/interpolation/element/Quad2D.h"
#include "atlas/util/Point.h"

#include "tests/AtlasTestEnvironment.h"

using atlas::PointXY;
using atlas::interpolation::element::Quad2D;
using atlas::interpolation::method::Intersect;

namespace atlas {
namespace test {

//----------------------------------------------------------------------------------------------------------------------

const double relative_error = 0.0001;

CASE("test_quad_area") {
    PointXY v0(0., 0.);
    PointXY v1(1., 0.);
    PointXY v2(1., 1.);
    PointXY v3(0., 1.);

    Quad2D quad1(v0.data(), v1.data(), v2.data(), v3.data());

    EXPECT(quad1.validate());

    double area = quad1.area();

    std::cout << "area " << area << std::endl;
    EXPECT_APPROX_EQ(area, 1.0, relative_error);

    PointXY c0(-2., -2.);  // 4
    PointXY c1(3., -2.);   // 6
    PointXY c2(3., 0.5);   // 1.5
    PointXY c3(-2., 0.5);  // 1

    Quad2D quad2(c0.data(), c1.data(), c2.data(), c3.data());

    EXPECT(quad2.validate());

    area = quad2.area();

    std::cout << "area " << area << std::endl;
    EXPECT_APPROX_EQ(area, 12.5, relative_error);
}

CASE("test_quadrilateral_intersection_refquad") {
    PointXY v0(0., 0.);
    PointXY v1(1., 0.);
    PointXY v2(1., 1.);
    PointXY v3(0., 1.);

    Quad2D quad(v0.data(), v1.data(), v2.data(), v3.data());

    EXPECT(quad.validate());

    PointLonLat orig(0.25, 0.25);

    Intersect isect = quad.intersects(orig);

    EXPECT(isect);
    EXPECT_APPROX_EQ(isect.u, 0.25, relative_error);
    EXPECT_APPROX_EQ(isect.v, 0.25, relative_error);
}

CASE("test_quadrilateral_remap_refquad") {
    PointXY v0(0., 0.);
    PointXY v1(1., 0.);
    PointXY v2(1., 1.);
    PointXY v3(0., 1.);

    Quad2D quad(v0.data(), v1.data(), v2.data(), v3.data());

    EXPECT(quad.validate());

    PointXY orig(0.25, 0.25);

    Intersect isect = quad.localRemap(orig);

    EXPECT(isect);
    EXPECT_APPROX_EQ(isect.u, 0.25, relative_error);
    EXPECT_APPROX_EQ(isect.v, 0.25, relative_error);
}

CASE("test_quadrilateral_intersection_doublequad") {
    PointXY v0(0., 0.);
    PointXY v1(2., 0.);
    PointXY v2(2., 2.);
    PointXY v3(0., 2.);

    Quad2D quad(v0.data(), v1.data(), v2.data(), v3.data());

    EXPECT(quad.validate());

    PointXY orig(0.5, 0.5);

    Intersect isect = quad.localRemap(orig);

    EXPECT(isect);
    EXPECT_APPROX_EQ(isect.u, 0.25, relative_error);
    EXPECT_APPROX_EQ(isect.v, 0.25, relative_error);
}

CASE("test_quadrilateral_intersection_rotatedquad") {
    PointXY v0(0., -1.);
    PointXY v1(1., 0.);
    PointXY v2(0., 1.);
    PointXY v3(-1., 0.);

    Quad2D quad(v0.data(), v1.data(), v2.data(), v3.data());

    EXPECT(quad.validate());

    PointXY orig(0., 0.);

    Intersect isect = quad.localRemap(orig);

    EXPECT(isect);
    EXPECT_APPROX_EQ(isect.u, 0.5, relative_error);
    EXPECT_APPROX_EQ(isect.v, 0.5, relative_error);
}

CASE("test_quadrilateral_intersection_arbitrary") {
    PointXY v0(338.14, 54.6923);
    PointXY v1(340.273, 54.6778);
    PointXY v2(340.312, 55.9707);
    PointXY v3(338.155, 55.9852);

    Quad2D quad(v0.data(), v1.data(), v2.data(), v3.data());

    EXPECT(quad.validate());

    PointXY orig(339, 55);

    Intersect isect = quad.localRemap(orig);

    std::cout << isect.u << " " << isect.v << std::endl;
    EXPECT(isect);
    EXPECT_APPROX_EQ(isect.u, 0.400390, relative_error);
    EXPECT_APPROX_EQ(isect.v, 0.242483, relative_error);
}

CASE("test_quadrilateral_intersection_nointersect") {
    PointXY v0(0., -1.);
    PointXY v1(1., 0.);
    PointXY v2(0., 1.);
    PointXY v3(-1., 0.);

    Quad2D quad(v0.data(), v1.data(), v2.data(), v3.data());

    EXPECT(quad.validate());

    PointXY orig(2., 2.);

    Intersect isect = quad.localRemap(orig);
    EXPECT(!isect);
}

CASE("test_quadrilateral_intersection_corners") {
    PointXY v0(0.0, -2.0);
    PointXY v1(2.5, 0.0);
    PointXY v2(0.0, 3.5);
    PointXY v3(-1.5, 0.0);

    Quad2D quad(v0.data(), v1.data(), v2.data(), v3.data());

    EXPECT(quad.validate());

    std::vector<PointXY> corners;
    corners.emplace_back(0.0, -2.0);
    corners.emplace_back(2.5, 0.0);
    corners.emplace_back(0.0, 3.5);
    corners.emplace_back(-1.5, 0.0);

    std::vector<std::pair<double, double>> uvs;
    uvs.emplace_back(0., 0.);
    uvs.emplace_back(1., 0.);
    uvs.emplace_back(1., 1.);
    uvs.emplace_back(0., 1.);

    for (size_t i = 0; i < 4; ++i) {
        PointXY orig = corners[i];

        Intersect isect = quad.localRemap(orig);

        EXPECT(isect);
        EXPECT_APPROX_EQ(isect.u, uvs[i].first, relative_error);
        EXPECT_APPROX_EQ(isect.v, uvs[i].second, relative_error);
    }
}

//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
