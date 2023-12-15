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
#include "atlas/util/Constants.h"
#include "atlas/util/function/VortexRollup.h"

#include "tests/AtlasTestEnvironment.h"

namespace atlas {
namespace test {

using namespace atlas::util;
using namespace atlas::array::helpers;

constexpr auto Rank2dField = 2;
constexpr auto Rank3dField = 3;

// Return (u, v) field with vortex_rollup as the streamfunction.
// This has no physical significance, but it makes a nice swirly field.
std::pair<double, double> vortexHorizontal(double lon, double lat) {

  // set hLon and hLat step size.
  const double hLon = 0.0001;
  const double hLat = 0.0001;

  // Get finite differences.

  // Set u.
  const double u = (function::vortex_rollup(lon, lat + 0.5 * hLat, 0.1) -
                    function::vortex_rollup(lon, lat - 0.5 * hLat, 0.1)) /
                   hLat;

  const double v = -(function::vortex_rollup(lon + 0.5 * hLon, lat, 0.1) -
                     function::vortex_rollup(lon - 0.5 * hLon, lat, 0.1)) /
                   (hLon * std::cos(lat * util::Constants::degreesToRadians()));

  return std::make_pair(u, v);
}

double vortexVertical(double lon, double lat) {
  return function::vortex_rollup(lon, lat, 0.1);
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

// Helper function to generate a NodeColumns functionspace
const auto generateNodeColums(const std::string& gridName,
                              const std::string& meshName) {
  const auto grid = Grid(gridName);
  const auto mesh = MeshGenerator(meshName).generate(grid);
  return functionspace::NodeColumns(mesh);
}

// Helper struct to key different Functionspaces to strings
struct FunctionSpaceFixtures {
  static const FunctionSpace& get(const std::string& fixture) {
    static std::map<std::string_view, FunctionSpace> functionSpaces = {
        {"cubedsphere_mesh",
         generateNodeColums("CS-LFR-48", "cubedsphere_dual")},
        {"gaussian_mesh", generateNodeColums("O48", "structured")}};
    return functionSpaces.at(fixture);
  }
};

// Helper struct to key different grid configs to strings
struct FieldSpecFixtures {
  static const Config& get(const std::string& fixture) {
    static std::map<std::string_view, Config> fieldSpecs = {
        {"2vector", option::name("test field") | option::variables(2) |
                        option::type("vector")},
        {"3vector", option::name("test field") | option::variables(3) |
                        option::type("vector")}};
    return fieldSpecs.at(fixture);
  }
};


template <int Rank>
void testInterpolation(const Config& config) {

  const auto& sourceFunctionSpace =
      FunctionSpaceFixtures::get(config.getString("source_fixture"));
  const auto& targetFunctionSpace =
      FunctionSpaceFixtures::get(config.getString("target_fixture"));

  auto sourceFieldSet = FieldSet{};
  auto targetFieldSet = FieldSet{};

  const auto sourceLonLat =
      array::make_view<double, 2>(sourceFunctionSpace.lonlat());
  const auto targetLonLat =
      array::make_view<double, 2>(targetFunctionSpace.lonlat());

  auto fieldSpec =
      FieldSpecFixtures::get(config.getString("field_spec_fixture"));
  if constexpr (Rank == 3) fieldSpec.set("levels", 2);

  auto sourceField = array::make_view<double, Rank>(
      sourceFieldSet.add(sourceFunctionSpace.createField<double>(fieldSpec)));

  auto targetField = array::make_view<double, Rank>(
      targetFieldSet.add(targetFunctionSpace.createField<double>(fieldSpec)));

  ArrayForEach<0>::apply(std::tie(sourceLonLat, sourceField),
                         [](auto&& lonLat, auto&& sourceColumn) {

    const auto setElems = [&](auto&& sourceElem) {
      std::tie(sourceElem(0), sourceElem(1)) =
          vortexHorizontal(lonLat(0), lonLat(1));
      if (sourceElem.size() == 3) {
        sourceElem(2) = vortexVertical(lonLat(0), lonLat(1));
      }
    };
    if constexpr (Rank == 2) { setElems(sourceColumn); }
    else if constexpr (Rank == 3) {
        ArrayForEach<0>::apply(std::tie(sourceColumn), setElems);
    }
  });
  sourceFieldSet.set_dirty(false);

  const auto interp = Interpolation(config.getSubConfiguration("scheme"),
                                    sourceFunctionSpace, targetFunctionSpace);

  interp.execute(sourceFieldSet, targetFieldSet);
  targetFieldSet.haloExchange();

  auto errorFieldSpec = fieldSpec;
  errorFieldSpec.remove("variables");

  auto errorField = array::make_view<double, Rank - 1>(targetFieldSet.add(
      targetFunctionSpace.createField<double>(errorFieldSpec)));

  auto maxError = 0.;
  ArrayForEach<0>::apply(std::tie(targetLonLat, targetField, errorField),
                         [&](auto&& lonLat, auto&& targetColumn,
                             auto&& errorColumn) {

    const auto calcError = [&](auto&& targetElem, auto&& errorElem) {
      auto trueValue = std::vector<double>(targetElem.size());
      std::tie(trueValue[0], trueValue[1]) =
          vortexHorizontal(lonLat(0), lonLat(1));
      if (targetElem.size() == 3) {
        trueValue[2] = vortexVertical(lonLat(0), lonLat(1));
      }

      auto errorSqrd = 0.;
      for (auto k = 0; k < targetElem.size(); ++k) {
        errorSqrd +=
            (targetElem(k) - trueValue[k]) * (targetElem(k) - trueValue[k]);
      }

      errorElem = std::sqrt(errorSqrd);
      maxError = std::max(maxError, static_cast<double>(errorElem));
    };

    if constexpr(Rank == 2) { calcError(targetColumn, errorColumn); }
    else if constexpr (Rank == 3) {
        ArrayForEach<0>::apply(std::tie(targetColumn, errorColumn), calcError);
    }
  });

  EXPECT_APPROX_EQ(maxError, 0., config.getDouble("tol"));

  gmshOutput(config.getString("file_id") + "_source.msh", sourceFieldSet);
  gmshOutput(config.getString("file_id") + "_target.msh", targetFieldSet);
}

CASE("cubed sphere vector interpolation (3d-field, 2-vector)") {

  const auto baseInterpScheme = util::Config("type", "cubedsphere-bilinear");
  const auto interpScheme =
      util::Config("type", "spherical-vector").set("scheme", baseInterpScheme);
  const auto cubedSphereConf = Config("source_fixture", "cubedsphere_mesh")
                                   .set("target_fixture", "gaussian_mesh")
                                   .set("field_spec_fixture", "2vector")
                                   .set("file_id", "spherical_vector_cs2")
                                   .set("scheme", interpScheme)
                                   .set("tol", 0.00018);

  testInterpolation<Rank3dField>((cubedSphereConf));
}

CASE("cubed sphere vector interpolation (3d-field, 3-vector)") {

  const auto baseInterpScheme = util::Config("type", "cubedsphere-bilinear");
  const auto interpScheme =
      util::Config("type", "spherical-vector").set("scheme", baseInterpScheme);
  const auto cubedSphereConf = Config("source_fixture", "cubedsphere_mesh")
                                   .set("target_fixture", "gaussian_mesh")
                                   .set("field_spec_fixture", "3vector")
                                   .set("file_id", "spherical_vector_cs3")
                                   .set("scheme", interpScheme)
                                   .set("tol", 0.00096);

  testInterpolation<Rank3dField>((cubedSphereConf));
}

CASE("finite element vector interpolation (2d-field, 2-vector)") {

  const auto baseInterpScheme = util::Config("type", "finite-element");
  const auto interpScheme =
      util::Config("type", "spherical-vector").set("scheme", baseInterpScheme);
  const auto cubedSphereConf = Config("source_fixture", "gaussian_mesh")
                                   .set("target_fixture", "cubedsphere_mesh")
                                   .set("field_spec_fixture", "2vector")
                                   .set("file_id", "spherical_vector_fe")
                                   .set("scheme", interpScheme)
                                   .set("tol", 0.00015);

  testInterpolation<Rank2dField>((cubedSphereConf));
}
}
}

int main(int argc, char** argv) { return atlas::test::run(argc, argv); }
