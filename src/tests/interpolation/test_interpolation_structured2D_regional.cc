/*
 * (C) Copyright 2024 Meteorologisk Institutt
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "atlas/array.h"
#include "atlas/field.h"
#include "atlas/functionspace/StructuredColumns.h"
#include "atlas/grid.h"
#include "atlas/interpolation.h"
#include "atlas/option.h"
#include "atlas/util/Config.h"

#include "tests/AtlasTestEnvironment.h"

using atlas::functionspace::StructuredColumns;
using atlas::util::Config;

namespace atlas {
namespace test {

//-----------------------------------------------------------------------------

static Config gridConfigs() {
    Config gridConfigs;

    Config projectionConfig;
    projectionConfig.set("type", "lambert_conformal_conic");
    projectionConfig.set("latitude0", 56.3);
    projectionConfig.set("longitude0", 0.0);

    Config gridConfig;
    gridConfig.set("type", "regional");
    std::vector<double> lonlat = {9.9, 56.3};
    gridConfig.set("lonlat(centre)", lonlat);
    gridConfig.set("projection", projectionConfig);

    const size_t sourceNx = 21;
    const size_t sourceNy = 31;
    const double sourceDx = 8.0e3;
    const double sourceDy = 9.0e3;

    const size_t xFactor = 4;
    const size_t yFactor = 3;

    const size_t targetNx = (sourceNx-1)*xFactor+1;
    const size_t targetNy = (sourceNy-1)*yFactor+1;
    const double targetDx = static_cast<double>(sourceNx-1)/static_cast<double>(targetNx-1)*sourceDx;
    const double targetDy = static_cast<double>(sourceNy-1)/static_cast<double>(targetNy-1)*sourceDy;

    gridConfig.set("nx", sourceNx);
    gridConfig.set("ny", sourceNy);
    gridConfig.set("dx", sourceDx);
    gridConfig.set("dy", sourceDy);
    gridConfigs.set("source", gridConfig);
    gridConfigs.set("source normalization", (sourceNx-1)*(sourceNy-1));

    gridConfig.set("nx", targetNx);
    gridConfig.set("ny", targetNy);
    gridConfig.set("dx", targetDx);
    gridConfig.set("dy", targetDy);
    gridConfigs.set("target", gridConfig);
    gridConfigs.set("target normalization", (targetNx-1)*(targetNy-1));

    return gridConfigs;
}

// Dot product
double dotProd(const Field& field01, const Field& field02) {
  double dprod{};

  const size_t ndim = field01.levels() > 0 ? 2 : 1;
  if (ndim == 1) {
    const auto field01_view = array::make_view<double, 1>(field01);
    const auto field02_view = array::make_view<double, 1>(field02);
    for (idx_t i=0; i < field01_view.shape(0); ++i) {
      dprod += field01_view(i) * field02_view(i);
    }
  } else if (ndim == 2) {
    const auto field01_view = array::make_view<double, 2>(field01);
    const auto field02_view = array::make_view<double, 2>(field02);
    for (idx_t l=0; l < field01_view.shape(1); ++l) {
      for (idx_t i=0; i < field01_view.shape(0); ++i) {
        dprod += field01_view(i, l) * field02_view(i, l);
      }
    }
  }
  eckit::mpi::comm().allReduceInPlace(dprod, eckit::mpi::Operation::SUM);

  return dprod;
}


CASE("test_interpolation_structured2D_regional_1d") {
    Grid sourceGrid(gridConfigs().getSubConfiguration("source"));
    Grid targetGrid(gridConfigs().getSubConfiguration("target"));

    StructuredColumns sourceFs(sourceGrid, option::halo(1));
    StructuredColumns targetFs(targetGrid, option::halo(1));

    Interpolation interpolation(Config("type", "regional-bilinear"), sourceFs, targetFs);

    auto sourceField = sourceFs.createField<double>(Config("name", "source"));
    auto targetField = targetFs.createField<double>(Config("name", "target"));

    // Accuracy test
    const auto sourceIView = array::make_view<idx_t, 1>(sourceFs.index_i());
    const auto sourceJView = array::make_view<idx_t, 1>(sourceFs.index_j());
    auto sourceView = array::make_view<double, 1>(sourceField);
    const auto sourceGhostView = atlas::array::make_view<int, 1>(sourceFs.ghost());
    sourceView.assign(0.0);
    for (idx_t i = 0; i < sourceFs.size(); ++i) {
      if (sourceGhostView(i) == 0) {
        sourceView(i) = static_cast<double>((sourceIView(i)-1)*(sourceJView(i)-1))
                       /static_cast<double>(gridConfigs().getInt("source normalization"));
      }
    }

    interpolation.execute(sourceField, targetField);

    const auto targetIView = array::make_view<idx_t, 1>(targetFs.index_i());
    const auto targetJView = array::make_view<idx_t, 1>(targetFs.index_j());
    const auto targetView = array::make_view<double, 1>(targetField);
    const auto targetGhostView = atlas::array::make_view<int, 1>(targetFs.ghost());
    const double tolerance = 1.e-12;
    for (idx_t i = 0; i < targetFs.size(); ++i) {
      if (targetGhostView(i) == 0) {
        const double targetTest = static_cast<double>((targetIView(i)-1)*(targetJView(i)-1))
                                 /static_cast<double>(gridConfigs().getInt("target normalization"));
        EXPECT_APPROX_EQ(targetView(i), targetTest, tolerance);
      }
    }

    // Adjoint test
    auto targetAdjoint = targetFs.createField<double>();
    array::make_view<double, 1>(targetAdjoint).assign(array::make_view<double, 1>(targetField));
    targetAdjoint.adjointHaloExchange();

    auto sourceAdjoint = sourceFs.createField<double>();
    array::make_view<double, 1>(sourceAdjoint).assign(0.);
    interpolation.execute_adjoint(sourceAdjoint, targetAdjoint);

    const auto yDotY    = dotProd(targetField, targetField);
    const auto xDotXAdj = dotProd(sourceField, sourceAdjoint);

    EXPECT_APPROX_EQ(yDotY / xDotXAdj, 1., 1e-14);
}


CASE("test_interpolation_structured2D_regional_2d") {
    Grid sourceGrid(gridConfigs().getSubConfiguration("source"));
    Grid targetGrid(gridConfigs().getSubConfiguration("target"));

    const idx_t nlevs = 2;
    StructuredColumns sourceFs(sourceGrid, option::halo(1) | option::levels(nlevs));
    StructuredColumns targetFs(targetGrid, option::halo(1) | option::levels(nlevs));

    Interpolation interpolation(Config("type", "regional-bilinear"), sourceFs, targetFs);

    auto sourceField = sourceFs.createField<double>(Config("name", "source"));
    auto targetField = targetFs.createField<double>(Config("name", "target"));

    // Accuracy test
    const auto sourceIView = array::make_view<idx_t, 1>(sourceFs.index_i());
    const auto sourceJView = array::make_view<idx_t, 1>(sourceFs.index_j());
    auto sourceView = array::make_view<double, 2>(sourceField);
    const auto sourceGhostView = atlas::array::make_view<int, 1>(sourceFs.ghost());
    sourceView.assign(0.0);
    for (idx_t i = 0; i < sourceFs.size(); ++i) {
      if (sourceGhostView(i) == 0) {
        for (idx_t k = 0; k < nlevs; ++k) {
          sourceView(i, k) = static_cast<double>((sourceIView(i)-1)*(sourceJView(i)-1))
                           /static_cast<double>(gridConfigs().getInt("source normalization"));
        }
      }
    }

    interpolation.execute(sourceField, targetField);

    const auto targetIView = array::make_view<idx_t, 1>(targetFs.index_i());
    const auto targetJView = array::make_view<idx_t, 1>(targetFs.index_j());
    const auto targetView = array::make_view<double, 2>(targetField);
    const auto targetGhostView = atlas::array::make_view<int, 1>(targetFs.ghost());
    const double tolerance = 1.e-12;
    for (idx_t i = 0; i < targetFs.size(); ++i) {
      if (targetGhostView(i) == 0) {
        const double targetTest = static_cast<double>((targetIView(i)-1)*(targetJView(i)-1))
                                 /static_cast<double>(gridConfigs().getInt("target normalization"));
        for (idx_t k = 0; k < nlevs; ++k) {
          EXPECT_APPROX_EQ(targetView(i, k), targetTest, tolerance);
        }
      }
    }

    // Adjoint test
    auto targetAdjoint = targetFs.createField<double>();
    array::make_view<double, 2>(targetAdjoint).assign(array::make_view<double, 2>(targetField));
    targetAdjoint.adjointHaloExchange();

    auto sourceAdjoint = sourceFs.createField<double>();
    array::make_view<double, 2>(sourceAdjoint).assign(0.);
    interpolation.execute_adjoint(sourceAdjoint, targetAdjoint);

    const auto yDotY    = dotProd(targetField, targetField);
    const auto xDotXAdj = dotProd(sourceField, sourceAdjoint);

    EXPECT_APPROX_EQ(yDotY / xDotXAdj, 1., 1e-14);
}

}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
