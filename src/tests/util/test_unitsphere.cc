/*
 * (C) Crown Copyright 2023 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "atlas/util/Point.h"
#include "atlas/util/Geometry.h"

#include "tests/AtlasTestEnvironment.h"

namespace atlas {
namespace test {

using namespace atlas::util;

CASE("great-circle course") {

  geometry::UnitSphere g;

  const auto point1 = PointLonLat(-71.6, -33.0);  // Valparaiso
  const auto point2 = PointLonLat(121.8, 31.4);   // Shanghai
  const auto point3 = PointLonLat(0., 89.);
  const auto point4 = PointLonLat(180., 89.);

  const auto targetCourse1 = -94.41;
  const auto targetCourse2 = -78.42;
  const auto targetCourse3 = 0.;
  const auto targetCourse4 = 180.;

  const auto[ course1, course2 ] = g.greatCircleCourse(point1, point2);
  const auto[ course3, course4 ] = g.greatCircleCourse(point3, point4);

  EXPECT_APPROX_EQ(course1, targetCourse1, 0.01);
  EXPECT_APPROX_EQ(course2, targetCourse2, 0.01);
  EXPECT_APPROX_EQ(course3, targetCourse3, 1.e-14);
  EXPECT_APPROX_EQ(std::abs(course4), targetCourse4, 1.e-14);
}
} // namespace test
} // namespace atlas

int main(int argc, char** argv) { return atlas::test::run(argc, argv); }
