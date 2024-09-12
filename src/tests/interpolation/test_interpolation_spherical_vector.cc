/*
 * (C) Crown Copyright 2024 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <cmath>
#include <limits>
#include <utility>

#include "atlas/array.h"
#include "atlas/array/helpers/ArrayForEach.h"
#include "atlas/field.h"
#include "atlas/functionspace.h"
#include "atlas/grid.h"
#include "atlas/interpolation.h"
#include "atlas/mesh.h"
#include "atlas/meshgenerator.h"
#include "atlas/output/Gmsh.h"
#include "atlas/util/Config.h"
#include "atlas/util/Constants.h"
#include "atlas/util/Point.h"
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
  const auto functionSpace = fieldSet[0].functionspace();
  const auto structuredColumns = functionspace::StructuredColumns(functionSpace);
  const auto nodeColumns = functionspace::NodeColumns(functionSpace);
  const auto mesh =
      structuredColumns ? Mesh(structuredColumns.grid()) : nodeColumns.mesh();

  const auto gmshConfig = Config("coordinates", "xyz") | Config("ghost", true) |
                          Config("info", true);
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

// Helper function to create part-empty PointCloud
const auto generateEmptyPointCloud() {
  const auto functionSpace = functionspace::PointCloud(std::vector<PointXY>{});
  return functionSpace;
}

// Helper struct to key different Functionspaces to strings
struct FunctionSpaceFixtures {
  static const FunctionSpace& get(const std::string& fixture) {
    static auto functionSpaces =
        std::map<std::string_view, FunctionSpace>{
            {"cubedsphere_mesh",
             generateNodeColums("CS-LFR-48", "cubedsphere_dual")},
            {"gaussian_mesh", generateNodeColums("O48", "structured")},
            {"structured_columns",
             functionspace::StructuredColumns(Grid("O48"), option::halo(1))},
            {"structured_columns_classic",
             functionspace::StructuredColumns(Grid("F48"), option::halo(1))},
            {"structured_columns_classic_halo2",
             functionspace::StructuredColumns(Grid("F48"), option::halo(2))},
            {"structured_columns_classic_highres_halo2",
             functionspace::StructuredColumns(Grid("F96"), option::halo(2))},
            {"structured_columns_halo2",
             functionspace::StructuredColumns(Grid("O48"), option::halo(2))},
            {"structured_columns_lowres",
             functionspace::StructuredColumns(Grid("O24"), option::halo(1))},
            {"structured_columns_hires",
             functionspace::StructuredColumns(Grid("O96"), option::halo(1))},
            {"empty_point_cloud", generateEmptyPointCloud()}};
    return functionSpaces.at(fixture);
  }
};

// Helper struct to key different grid configs to strings
struct FieldSpecFixtures {
  static const Config& get(const std::string& fixture) {
    static const auto fieldSpecs = std::map<std::string_view, Config>{
        {"2vector", option::name("test field") | option::variables(2) |
                        option::type("vector")},
        {"3vector", option::name("test field") | option::variables(3) |
                        option::type("vector")}};
    return fieldSpecs.at(fixture);
  }
};

// Helper struct to key different interpolation schemes to strings
struct InterpSchemeFixtures {
  static const Config& get(const std::string& fixture) {
    static const auto cubedsphereBilinear =
        option::type("cubedsphere-bilinear") | Config("adjoint", true);
    static const auto finiteElement =
        option::type("finite-element") | Config("adjoint", true);
    static const auto structuredLinear = option::type("structured-linear2D") |
                                         option::halo(1) |
                                         Config("adjoint", true);
    static const auto structuredCubic = option::type("structured-bicubic") |
                                        option::halo(2) |
                                        Config("adjoint", true);
    static const auto sphericalVector =
        option::type("spherical-vector") | Config("adjoint", true);

    static const auto interpSchemes = std::map<std::string_view, Config>{
        {"cubedsphere_bilinear", cubedsphereBilinear},
        {"finite_element", finiteElement},
        {"structured_linear", structuredLinear},
        {"structured_cubic", structuredCubic},
        {"cubedsphere_bilinear_spherical",
         sphericalVector | Config("scheme", cubedsphereBilinear)},
        {"finite_element_spherical",
         sphericalVector | Config("scheme", finiteElement)},
        {"structured_linear_spherical",
         sphericalVector | Config("scheme", structuredLinear)},
        {"structured_cubic_spherical",
         sphericalVector | Config("scheme", structuredCubic)}};
    return interpSchemes.at(fixture);
  }
};

template <int Rank>
double dotProduct(const array::ArrayView<double, Rank>& a,
                  const array::ArrayView<double, Rank>& b) {
  auto dotProd = 0.;
  arrayForEachDim(std::make_integer_sequence<int, Rank>{}, std::tie(a, b),
                  [&](const double& aElem, const double& bElem) {
                    dotProd += aElem * bElem;
                  });
  return dotProd;
}

template <int Rank>
int countNans(const array::ArrayView<double, Rank>& view) {
  auto nNans = 0;
  arrayForEachDim(std::make_integer_sequence<int, Rank>{}, std::tie(view),
                  [&](const double& viewElem) {
                    if (!std::isfinite(viewElem)) {
                      ++nNans;
                    }
                  });
  return nNans;
}

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
  if constexpr (Rank == 3) {
      fieldSpec.set("levels", 2);
  }

  auto sourceField =
      sourceFieldSet.add(sourceFunctionSpace.createField<double>(fieldSpec));
  auto targetField =
      targetFieldSet.add(targetFunctionSpace.createField<double>(fieldSpec));

  auto sourceView = array::make_view<double, Rank>(sourceField);
  auto targetView = array::make_view<double, Rank>(targetField);

  ArrayForEach<0>::apply(
      std::tie(sourceLonLat, sourceView),
      [](auto&& lonLat, auto&& sourceColumn) {
        const auto setElems = [&](auto&& sourceElem) {
          std::tie(sourceElem(0), sourceElem(1)) =
              vortexHorizontal(lonLat(0), lonLat(1));
          if (sourceElem.size() == 3) {
            sourceElem(2) = vortexVertical(lonLat(0), lonLat(1));
          }
        };
        if constexpr (Rank == 2) {
          setElems(sourceColumn);
        } else if constexpr (Rank == 3) {
          ArrayForEach<0>::apply(std::tie(sourceColumn), setElems);
        }
      });

  const auto interp = Interpolation(
      InterpSchemeFixtures::get(config.getString("interp_fixture")),
      sourceFunctionSpace, targetFunctionSpace);

  interp.execute(sourceFieldSet, targetFieldSet);
  targetFieldSet.haloExchange();

  auto errorFieldSpec = fieldSpec;
  errorFieldSpec.remove("variables").set("name", "error field");

  auto errorView = array::make_view<double, Rank - 1>(targetFieldSet.add(
      targetFunctionSpace.createField<double>(errorFieldSpec)));
  errorView.assign(0.);

  if (config.has("tol")) {
    auto maxError = 0.;
    ArrayForEach<0>::apply(
        std::tie(targetLonLat, targetView, errorView),
        [&](auto&& lonLat, auto&& targetColumn, auto&& errorColumn) {
          const auto calcError = [&](auto&& targetElem, auto&& errorElem) {
            auto trueValue = std::vector<double>(targetElem.size());
            std::tie(trueValue[0], trueValue[1]) =
                vortexHorizontal(lonLat(0), lonLat(1));
            if (targetElem.size() == 3) {
              trueValue[2] = vortexVertical(lonLat(0), lonLat(1));
            }

            auto errorSqrd = 0.;
            for (auto k = 0; k < targetElem.size(); ++k) {
              errorSqrd += (targetElem(k) - trueValue[k]) *
                           (targetElem(k) - trueValue[k]);
            }

            errorElem = std::sqrt(errorSqrd);
            maxError = std::max(maxError, static_cast<double>(errorElem));
          };

          if constexpr (Rank == 2) {
            calcError(targetColumn, errorColumn);
          }
          else if constexpr (Rank == 3) {
            ArrayForEach<0>::apply(std::tie(targetColumn, errorColumn),
                                   calcError);
          }
        });

    EXPECT_APPROX_EQ(maxError, 0., config.getDouble("tol"));
  }

  if (config.has("file_id")) {
    gmshOutput(config.getString("file_id") + "_source.msh", sourceFieldSet);
    gmshOutput(config.getString("file_id") + "_target.msh", targetFieldSet);
  }

  // Adjoint test
  auto targetAdjoint = targetFunctionSpace.createField<double>(fieldSpec);
  auto targetAdjointView = array::make_view<double, Rank>(targetAdjoint);
  targetAdjoint.array().copy(targetField);
  targetAdjoint.adjointHaloExchange();

  auto sourceAdjoint = sourceFunctionSpace.createField<double>(fieldSpec);
  auto sourceAdjointView = array::make_view<double, Rank>(sourceAdjoint);
  sourceAdjointView.assign(0.);

  sourceAdjoint.set_dirty(false);
  interp.execute_adjoint(sourceAdjoint, targetAdjoint);

  // Check fields for nans or +/-inf
  EXPECT_EQ(countNans(targetView), 0);
  EXPECT_EQ(countNans(sourceView), 0);
  EXPECT_EQ(countNans(targetAdjointView), 0);
  EXPECT_EQ(countNans(sourceAdjointView), 0);

  constexpr auto tinyNum = 1e-13;
  const auto targetDotTarget = dotProduct(targetView, targetView);
  const auto sourceDotSourceAdjoint = dotProduct(sourceView, sourceAdjointView);

  if (targetFunctionSpace.size() > 0) {
    const auto dotProdRatio = targetDotTarget / sourceDotSourceAdjoint;
    EXPECT_APPROX_EQ(dotProdRatio, 1., tinyNum);
  }
}


CASE("cubed sphere CS-LFR-48 vector interpolation (3d-field, 2-vector)") {
  const auto config =
      Config("source_fixture", "cubedsphere_mesh")
          .set("target_fixture", "gaussian_mesh")
          .set("field_spec_fixture", "2vector")
          .set("interp_fixture", "cubedsphere_bilinear_spherical")
          .set("file_id", "spherical_vector_cs2")
          .set("tol", 0.00018);

  testInterpolation<Rank3dField>((config));
}

CASE("cubed sphere CS-LFR-48 vector interpolation (3d-field, 3-vector)") {
  const auto config =
      Config("source_fixture", "cubedsphere_mesh")
          .set("target_fixture", "gaussian_mesh")
          .set("field_spec_fixture", "3vector")
          .set("interp_fixture", "cubedsphere_bilinear_spherical")
          .set("file_id", "spherical_vector_cs3")
          .set("tol", 0.00096);

  testInterpolation<Rank3dField>((config));
}

CASE("cubed sphere CS-LFR-48 (spherical vector) to empty point cloud") {
    const auto config =
        Config("source_fixture", "cubedsphere_mesh")
            .set("target_fixture", "empty_point_cloud")
            .set("field_spec_fixture", "2vector")
            .set("interp_fixture", "cubedsphere_bilinear_spherical");

    testInterpolation<Rank2dField>((config));
}

CASE("cubed sphere CS-LFR-48 to empty point cloud") {
  const auto config =
      Config("source_fixture", "cubedsphere_mesh")
          .set("target_fixture", "empty_point_cloud")
          .set("field_spec_fixture", "2vector")
          .set("interp_fixture", "cubedsphere_bilinear");

  testInterpolation<Rank2dField>((config));
}


CASE("finite element to empty point cloud") {
  const auto config = Config("source_fixture", "gaussian_mesh")
                          .set("target_fixture", "cubedsphere_mesh")
                          .set("field_spec_fixture", "2vector")
                          .set("interp_fixture", "finite_element");

  testInterpolation<Rank2dField>((config));
}

CASE("finite element vector interpolation (2d-field, 2-vector)") {
  const auto config = Config("source_fixture", "gaussian_mesh")
                          .set("target_fixture", "cubedsphere_mesh")
                          .set("field_spec_fixture", "2vector")
                          .set("interp_fixture", "finite_element_spherical")
                          .set("file_id", "spherical_vector_fe")
                          .set("tol", 0.00015);

  testInterpolation<Rank2dField>((config));
}

CASE("structured columns F48 cubic vector spherical interpolation (3d-field, 2-vector)") {
  const auto config =
      Config("source_fixture", "structured_columns_classic_halo2")
          .set("target_fixture", "cubedsphere_mesh")
          .set("field_spec_fixture", "2vector")
          .set("interp_fixture", "structured_cubic_spherical")
          .set("file_id", "spherical_cubic_vector_classic_sc")
          .set("tol", 0.0000082);

  testInterpolation<Rank3dField>((config));
}

CASE("structured columns F96 cubic vector spherical interpolation (2d-field, 2-vector)") {
  const auto config =
      Config("source_fixture", "structured_columns_classic_highres_halo2")
          .set("target_fixture", "cubedsphere_mesh")
          .set("field_spec_fixture", "2vector")
          .set("interp_fixture", "structured_cubic_spherical")
          .set("file_id", "spherical_2D_cubic_vector_highres_classic_sc")
          .set("tol", 0.000000425);

  testInterpolation<Rank2dField>((config));
}

CASE("structured columns F96 cubic vector spherical interpolation (3d-field, 2-vector)") {
  const auto config =
      Config("source_fixture", "structured_columns_classic_highres_halo2")
          .set("target_fixture", "cubedsphere_mesh")
          .set("field_spec_fixture", "2vector")
          .set("interp_fixture", "structured_cubic_spherical")
          .set("file_id", "spherical_3D_cubic_vector_highres_classic_sc")
          .set("tol", 0.00000085);

  testInterpolation<Rank3dField>((config));
}

CASE("structured columns O24 linear vector interpolation (2d-field, 2-vector)") {
  const auto config = Config("source_fixture", "structured_columns_lowres")
                          .set("target_fixture", "gaussian_mesh")
                          .set("field_spec_fixture", "2vector")
                          .set("interp_fixture", "structured_linear_spherical")
                          .set("file_id", "spherical_vector_sc_lr")
                          .set("tol", 0.00056);

  testInterpolation<Rank2dField>((config));
}

CASE("structured columns O48 cubic vector spherical interpolation (3d-field, 2-vector)") {
  const auto config =
      Config("source_fixture", "structured_columns_halo2")
          .set("target_fixture", "cubedsphere_mesh")
          .set("field_spec_fixture", "2vector")
          .set("interp_fixture", "structured_cubic_spherical")
          .set("file_id", "spherical_cubic_vector_sc3")
          .set("tol",  0.000007);

  testInterpolation<Rank3dField>((config));
}

CASE("structured columns O48 linear vector interpolation (2d-field, 2-vector)") {
  const auto config = Config("source_fixture", "structured_columns")
                          .set("target_fixture", "cubedsphere_mesh")
                          .set("field_spec_fixture", "2vector")
                          .set("interp_fixture", "structured_linear_spherical")
                          .set("file_id", "spherical_vector_sc")
                          .set("tol", 0.00017);

  testInterpolation<Rank2dField>((config));
}

CASE("structured columns O48 to empty point cloud") {
  const auto config = Config("source_fixture", "structured_columns")
                          .set("target_fixture", "empty_point_cloud")
                          .set("field_spec_fixture", "2vector")
                          .set("interp_fixture", "structured_linear");

  testInterpolation<Rank2dField>((config));
}

CASE("structured columns O96 vector interpolation (2d-field, 2-vector, hi-res)") {
  const auto config = Config("source_fixture", "structured_columns_hires")
                          .set("target_fixture", "gaussian_mesh")
                          .set("field_spec_fixture", "2vector")
                          .set("interp_fixture", "structured_linear_spherical")
                          .set("file_id", "spherical_vector_sc_hr")
                          .set("tol", 0.000044);

  testInterpolation<Rank2dField>((config));
}


/*
CASE("gauss linear vector interpolation (2d-field, 2-vector)") {
  const auto config =
      Config("source_fixture", "structured_columns_classic")
          .set("target_fixture", "cubedsphere_mesh")
          .set("field_spec_fixture", "2vector")
          .set("interp_fixture", "structured_linear")
          .set("file_id", "linear_vector_rank2_sc1")
          .set("tol", 0.00087);

  testInterpolation<Rank2dField>((config));
}

CASE("gauss linear vector interpolation (3d-field, 2-vector)") {
  const auto config =
      Config("source_fixture", "structured_columns_classic")
          .set("target_fixture", "cubedsphere_mesh")
          .set("field_spec_fixture", "2vector")
          .set("interp_fixture", "structured_linear")
          .set("file_id", "linear_vector_rank3_sc1")
          .set("tol", 0.00087);

  testInterpolation<Rank3dField>((config));
}
*/


/*
CASE("gauss cubic vector interpolation (3d-field, 2-vector)") {
  const auto config =
      Config("source_fixture", "structured_columns_halo2")
          .set("target_fixture", "cubedsphere_mesh")
          .set("field_spec_fixture", "2vector")
          .set("interp_fixture", "structured_cubic")
          .set("file_id", "cubic_vector_sc2")
          .set("tol", 0.000575);

  testInterpolation<Rank3dField>((config));
}
*/

}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) { return atlas::test::run(argc, argv); }
