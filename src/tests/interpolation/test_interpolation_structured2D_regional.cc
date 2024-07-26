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

#include "eckit/utils/MD5.h"

#include "tests/AtlasTestEnvironment.h"

using atlas::functionspace::StructuredColumns;
using atlas::util::Config;

namespace atlas {
namespace test {

//-----------------------------------------------------------------------------

CASE("test_interpolation_structured2D_regional") {
    Config projectionConfig;
    projectionConfig.set("type", "lambert_conformal_conic");
    projectionConfig.set("latitude0", 56.3);
    projectionConfig.set("longitude0", 0.0);

    Config gridConfig;
    gridConfig.set("type", "regional");
    std::vector<double> lonlat = {9.9, 56.3};
    gridConfig.set("lonlat(centre)", lonlat);
    gridConfig.set("projection", projectionConfig);

    gridConfig.set("nx", 21);
    gridConfig.set("ny", 31);
    gridConfig.set("dx", 8.0e3);
    gridConfig.set("dy", 9.0e3);
    Grid sourceGrid(gridConfig);

    gridConfig.set("nx", 81);
    gridConfig.set("ny", 91);
    gridConfig.set("dx", 2.0e3);
    gridConfig.set("dy", 3.0e3);
    Grid targetGrid(gridConfig);

    StructuredColumns sourceFs(sourceGrid, option::halo(1));
    StructuredColumns targetFs(targetGrid, option::halo(1));

    Interpolation interpolation(Config("type", "regional-linear-2d"), sourceFs, targetFs);

    auto sourceField = sourceFs.createField<double>(Config("name", "source"));
    auto targetField = targetFs.createField<double>(Config("name", "target"));

    const auto indexIView = array::make_view<int, 1>(sourceFs.index_i());
    const auto indexJView = array::make_view<int, 1>(sourceFs.index_j());
    auto sourceView = array::make_view<double, 1>(sourceField);
    for (idx_t i = 0; i < sourceFs.size(); ++i) {
      sourceView(i) = indexIView(i)*31.0+indexJView(i);
    }

    interpolation.execute(sourceField, targetField);

    eckit::MD5 hash;
    const auto targetView = array::make_view<double, 1>(targetField);
    for (idx_t i = 0; i < targetFs.size(); ++i) {
      hash.add(targetView(i));
    }
    if (eckit::mpi::comm().rank() == 0) {
      ASSERT(hash.digest() == "f7ca2ef899f96e390a13db614d9fa231");
    } else if (eckit::mpi::comm().rank() == 1) {
      ASSERT(hash.digest() == "77e819709fe711d93311d8359de9724b");
    }
}

}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
