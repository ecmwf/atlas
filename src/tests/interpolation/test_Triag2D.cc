/*
 * (C) Crown Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "eckit/types/FloatCompare.h"

#include "atlas/interpolation/element/Triag2D.h"
#include "atlas/util/Point.h"

#include "tests/AtlasTestEnvironment.h"

using atlas::PointXY;
using atlas::interpolation::element::Triag2D;
using atlas::interpolation::method::Intersect;

namespace atlas {
namespace test {

//----------------------------------------------------------------------------------------------------------------------

const double relative_error = 0.0001;

CASE("test_triag_area") {
    PointXY v0(0., 0.);
    PointXY v1(1., 0.);
    PointXY v2(1., 1.);

    Triag2D triangle1(v0.data(), v1.data(), v2.data());

    EXPECT(triangle1.validate());

    double area = triangle1.area();

    std::cout << "area " << area << std::endl;
    EXPECT_APPROX_EQ(area, 0.5, relative_error);

    PointXY c0(-2., -2.);  // 4
    PointXY c1(3., -2.);   // 6
    PointXY c2(3., 0.5);   // 1.5

    Triag2D triangle2(c0.data(), c1.data(), c2.data());

    EXPECT(triangle2.validate());

    area = triangle2.area();

    std::cout << "area " << area << std::endl;
    EXPECT_APPROX_EQ(area, 6.25, relative_error);
}

CASE("test_intersection_equilateral_triangle") {
    PointXY v0(0., 0.);
    PointXY v1(1., 0.);
    PointXY v2(0.5, std::sqrt(3.0) / 2.0);

    Triag2D triangle(v0.data(), v1.data(), v2.data());

    EXPECT(triangle.validate());

    PointXY orig(0.5, std::sqrt(3.0) / 6.0);

    Intersect isect = triangle.intersects(orig);

    EXPECT(isect);
    std::cout << "isect.u " << isect.u << std::endl;
    std::cout << "isect.v " << isect.v << std::endl;
    EXPECT_APPROX_EQ(isect.u, 1. / 3., relative_error);
    EXPECT_APPROX_EQ(isect.v, 1. / 3., relative_error);
}

CASE("test_intersection_right_angled_triangle") {
    PointXY v0(0., 0.);
    PointXY v1(3., 0.);
    PointXY v2(3., 4.);

    Triag2D triangle(v0.data(), v1.data(), v2.data());

    EXPECT(triangle.validate());

    // average of the coordinates is the centre
    PointXY orig(2.0, (4.0 / 3.0));

    Intersect isect = triangle.intersects(orig);

    EXPECT(isect);
    std::cout << "isect.u " << isect.u << std::endl;
    std::cout << "isect.v " << isect.v << std::endl;
    EXPECT_APPROX_EQ(isect.u, 1. / 3., relative_error);
    EXPECT_APPROX_EQ(isect.v, 1. / 3., relative_error);
}

CASE("test_intersection_offset_right_angled_triangle") {
    PointXY v0(0., 0.);
    PointXY v1(3., 0.);
    PointXY v2(3., 4.);

    Triag2D triangle(v0.data(), v1.data(), v2.data());

    EXPECT(triangle.validate());

    PointXY orig(1.5, 1);

    Intersect isect = triangle.intersects(orig);

    EXPECT(isect);
    std::cout << "isect.u " << isect.u << std::endl;
    std::cout << "isect.v " << isect.v << std::endl;
    EXPECT_APPROX_EQ(isect.u, 0.25, relative_error);
    EXPECT_APPROX_EQ(isect.v, 0.25, relative_error);
}

CASE("test_intersection_rotatedtriangle") {
    const double root2 = std::sqrt(2.);
    PointXY v0(0., 0.);
    PointXY v1(3. / root2, -3. / root2);
    PointXY v2(7. / root2, 1. / root2);

    Triag2D triangle(v0.data(), v1.data(), v2.data());

    EXPECT(triangle.validate());

    PointXY orig(10. / (3. * root2), -2. / (3. * root2));

    Intersect isect = triangle.intersects(orig);

    EXPECT(isect);
    std::cout << "isect.u " << isect.u << std::endl;
    std::cout << "isect.v " << isect.v << std::endl;
    EXPECT_APPROX_EQ(isect.u, 1. / 3., relative_error);
    EXPECT_APPROX_EQ(isect.v, 1. / 3., relative_error);
}

CASE("test_intersection_nointersect") {
    PointXY v0(0., -1.);
    PointXY v1(1., 0.);
    PointXY v2(0., 1.);

    Triag2D triangle(v0.data(), v1.data(), v2.data());

    EXPECT(triangle.validate());

    PointXY orig(2., 2.);

    Intersect isect = triangle.intersects(orig);
    EXPECT(!isect);
}

CASE("test_intersection_corners") {
    PointXY v0(0.0, -2.0);
    PointXY v1(2.5, 0.0);
    PointXY v2(0.0, 3.5);

    Triag2D triangle(v0.data(), v1.data(), v2.data());

    EXPECT(triangle.validate());

    std::vector<PointXY> corners;
    corners.emplace_back(0.0, -2.0);
    corners.emplace_back(2.5, 0.0);
    corners.emplace_back(0.0, 3.5);

    std::vector<std::pair<double, double>> uvs;
    uvs.emplace_back(0., 0.);
    uvs.emplace_back(1., 0.);
    uvs.emplace_back(0., 1.);

    for (size_t i = 0; i < 3; ++i) {
        PointXY orig = corners[i];

        Intersect isect = triangle.intersects(orig);

        EXPECT(isect);
        std::cout << "isect.u " << isect.u << std::endl;
        std::cout << "isect.v " << isect.v << std::endl;
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
