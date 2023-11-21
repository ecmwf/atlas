/*
 * (C) Crown Copyright 2023 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "atlas/array.h"
#include "atlas/array/helpers/ArrayForEach.h"
#include "atlas/field.h"
#include "atlas/functionspace.h"
#include "atlas/grid.h"
#include "atlas/interpolation.h"
#include "atlas/mesh.h"
#include "atlas/meshgenerator.h"
#include "atlas/output/Gmsh.h"
#include "atlas/util/Point.h"
#include "atlas/util/Config.h"
#include "atlas/util/function/VortexRollup.h"
#include "atlas/util/UnitSphere.h"

#include "tests/AtlasTestEnvironment.h"

namespace atlas {
namespace test {

using namespace atlas::util;
using namespace atlas::array::helpers;

// Return (u, v) field with vortex_rollup as the streamfunction.
// This has no physical significance, but it makes a nice swirly field.
std::pair<double, double> vortexField(double lon, double lat) {

  // set hLon and hLat step size.
  const double hLon = 0.0001;
  const double hLat = 0.0001;

  // Get finite differences.

  // Set u.
  const double u = (util::function::vortex_rollup(lon, lat + 0.5 * hLat, 0.1) -
                    util::function::vortex_rollup(lon, lat - 0.5 * hLat, 0.1)) /
                   hLat;

  const double v = -(util::function::vortex_rollup(lon + 0.5 * hLon, lat, 0.1) -
                     util::function::vortex_rollup(lon - 0.5 * hLon, lat, 0.1)) /
                   (hLon * std::cos(lat * util::Constants::degreesToRadians()));

  return std::make_pair(u, v);
}

void gmshOutput(const std::string& fileName, const FieldSet& fieldSet) {

  const auto& functionSpace = fieldSet[0].functionspace();
  const auto& mesh = functionspace::NodeColumns(functionSpace).mesh();

  const auto gmshConfig = util::Config("coordinates", "xyz") |
                          util::Config("ghost", true) |
                          util::Config("info", true);
  const auto gmsh = output::Gmsh(fileName, gmshConfig);
  gmsh.write(mesh);
  gmsh.write(fieldSet, functionSpace);
}

void testInterpolation(const Config& config) {

  const auto sourceGrid = Grid(config.getString("source_grid"));
  const auto sourceMesh =
      MeshGenerator(config.getString("source_mesh")).generate(sourceGrid);
  const auto sourceFunctionSpace = functionspace::NodeColumns(sourceMesh);

  const auto targetGrid = Grid(config.getString("target_grid"));
  const auto targetMesh =
      MeshGenerator(config.getString("target_mesh")).generate(targetGrid);
  const auto targetFunctionSpace = functionspace::NodeColumns(targetMesh);

  auto sourceFieldSet = FieldSet{};
  auto targetFieldSet = FieldSet{};

  const auto sourceLonLat =
      array::make_view<double, 2>(sourceFunctionSpace.lonlat());
  const auto targetLonLat =
      array::make_view<double, 2>(targetFunctionSpace.lonlat());

  auto sourceField = array::make_view<double, 3>(sourceFieldSet.add(
      sourceFunctionSpace.createField<double>(option::name("test field") |
                                              option::levels(1) |
                                              option::variables(2))));

  auto targetField = array::make_view<double, 3>(targetFieldSet.add(
      targetFunctionSpace.createField<double>(option::name("test field") |
                                              option::levels(1) |
                                              option::variables(2))));

  ArrayForEach<0>::apply(std::tie(sourceLonLat, sourceField),
                         [](auto&& lonLat, auto&& sourceColumn) {
    ArrayForEach<0>::apply(std::tie(sourceColumn), [&](auto&& sourceElem) {
      std::tie(sourceElem(0), sourceElem(1)) =
          vortexField(lonLat(0), lonLat(1));
    });
  });

  const auto interp = Interpolation(config.getSubConfiguration("scheme"),
                                    sourceFunctionSpace, targetFunctionSpace);

  interp.execute(sourceFieldSet, targetFieldSet);
  targetFieldSet.haloExchange();

  auto errorField = array::make_view<double, 2>(
      targetFieldSet.add(targetFunctionSpace.createField<double>(
          option::name("error field") | option::levels(1))));

  ArrayForEach<0>::apply(
      std::tie(targetLonLat, targetField, errorField),
      [](auto&& lonLat, auto&& targetColumn, auto&& errorColumn) {
        ArrayForEach<0>::apply(std::tie(targetColumn, errorColumn),
                               [&](auto&& targetElem, auto&& errorElem) {

          auto deltaVal = vortexField(lonLat(0), lonLat(1));
          deltaVal.first -= targetElem(0);
          deltaVal.second -= targetElem(1);

          errorElem = std::sqrt(deltaVal.first * deltaVal.first +
                                deltaVal.second * deltaVal.second);
        });
      });

  gmshOutput(config.getString("file_id") + "_source.msh", sourceFieldSet);
  gmshOutput(config.getString("file_id") + "_target.msh", targetFieldSet);
}

CASE("great-circle course") {

  const auto point1 = PointLonLat(-71.6, -33.0);  // Valparaiso
  const auto point2 = PointLonLat(121.8, 31.4);   // Shanghai
  const auto point3 = PointLonLat(0., 89.);
  const auto point4 = PointLonLat(180., 89.);

  const auto targetCourse1 = -94.41;
  const auto targetCourse2 = -78.42;
  const auto targetCourse3 = 0.;
  const auto targetCourse4 = 180.;

  const auto[ course1, course2 ] = greatCircleCourse(point1, point2);
  const auto[ course3, course4 ] = greatCircleCourse(point3, point4);

  EXPECT_APPROX_EQ(course1, targetCourse1, 0.01);
  EXPECT_APPROX_EQ(course2, targetCourse2, 0.01);
  EXPECT_APPROX_EQ(course3, targetCourse3, 1.e-7);
  EXPECT_APPROX_EQ(std::abs(course4), targetCourse4, 1.e-7);
}

CASE("cubed sphere vector interpolation") {

  const auto baseInterpScheme =
      util::Config("type", "cubedsphere-bilinear").set("adjoint", true);
  const auto interpScheme = util::Config("type", "parallel-transport")
                                .set("scheme", baseInterpScheme);
  const auto cubedSphereConf = Config("source_grid", "CS-LFR-48")
                                   .set("source_mesh", "cubedsphere_dual")
                                   .set("target_grid", "O48")
                                   .set("target_mesh", "structured")
                                   .set("file_id", "parallel_transport_cs")
                                   .set("scheme", interpScheme);

  testInterpolation((cubedSphereConf));
}

CASE("finite element vector interpolation") {

  const auto baseInterpScheme =
      util::Config("type", "finite-element").set("adjoint", true);
  const auto interpScheme = util::Config("type", "parallel-transport")
                                .set("scheme", baseInterpScheme);
  const auto cubedSphereConf = Config("source_grid", "O48")
                                   .set("source_mesh", "structured")
                                   .set("target_grid", "CS-LFR-48")
                                   .set("target_mesh", "cubedsphere_dual")
                                   .set("file_id", "parallel_transport_fe")
                                   .set("scheme", interpScheme);

  testInterpolation((cubedSphereConf));
}

}
}

int main(int argc, char** argv) { return atlas::test::run(argc, argv); }
