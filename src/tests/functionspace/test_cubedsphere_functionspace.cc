/*
 * (C) Crown Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "atlas/array.h"
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/functionspace/CellColumns.h"
#include "atlas/functionspace/CubedSphereColumns.h"
#include "atlas/field/FieldSet.h"
#include "atlas/grid.h"
#include "atlas/grid/CubedSphereGrid.h"
#include "atlas/grid/Tiles.h"
#include "atlas/mesh.h"
#include "atlas/meshgenerator.h"
#include "atlas/meshgenerator/detail/cubedsphere/CubedSphereUtility.h"
#include "atlas/option.h"
#include "atlas/util/CoordinateEnums.h"

#include "tests/AtlasTestEnvironment.h"

namespace atlas {
namespace test {

  double testFunction(double lon, double lat) {
    return std::sin(3 * lon * M_PI / 180) * std::sin(2 * lat * M_PI / 180);
  }

  template<typename BaseFunctionSpace>
  void testFunctionSpace(const functionspace::CubedSphereColumns<BaseFunctionSpace>& functionspace) {

    // Make field.
    auto field = functionspace.template createField<double>(
      util::Config("name", "test field"));
    auto fieldView = array::make_view<double, 1>(field);

    // Get view of lonlat.
    const auto lonLatView = array::make_view<double, 2>(functionspace.lonlat());

    // Get view of ghost field if NodeColumns, halo field if CellColumns.
    const auto ghostView = functionspace.type() == "NodeColumns" ?
        array::make_view<idx_t, 1>( functionspace.mesh().nodes().ghost() ) :
        array::make_view<idx_t, 1>( functionspace.mesh().nodes().halo() );

    // Loop over all non halo elements of test field.
    idx_t testFuncCallCount = 0;
    functionspace.for_each(
      [&](idx_t index, idx_t t, idx_t i, idx_t j) {

        // Make sure index matches ijt.
        EXPECT(index == functionspace.index(t, i, j));

        // Check that indices of "+" stencil are valid.
        const auto badIdx = functionspace.invalid_index();
        EXPECT(functionspace.index(t, i - 1, j    ) != badIdx);
        EXPECT(functionspace.index(t, i + 1, j    ) != badIdx);
        EXPECT(functionspace.index(t, i    , j - 1) != badIdx);
        EXPECT(functionspace.index(t, i    , j + 1) != badIdx);

        // Make sure we're avoiding halos.
        EXPECT(!ghostView(index));

        // Set field values.
        fieldView(index) = testFunction(lonLatView(index, LON), lonLatView(index, LAT));
        ++testFuncCallCount;

    });

    // Make sure call count is less than functionspace.size() as we skipped halos.
    EXPECT(testFuncCallCount < functionspace.size());

    // Perform halo exchange.
    functionspace.haloExchange(field);

    // Loop over elements including halo
    testFuncCallCount = 0;
    functionspace.for_each(
      [&](idx_t index, idx_t t, idx_t i, idx_t j) {

        // Make sure index matches ijt.
        EXPECT(index == functionspace.index(t, i, j));

        // Set field values.
        EXPECT(is_approximately_equal(
          fieldView(index), testFunction(lonLatView(index, LON), lonLatView(index, LAT))));
        ++testFuncCallCount;

    }, true);

    // Make sure call count is equal to functionspace.size().
    EXPECT(testFuncCallCount == functionspace.size());


  }

CASE("cubedsphere_mesh_functionspace") {

  // Set grid.
  const auto grid = Grid("CS-LFR-C-12");

  // Set mesh config.
  const auto meshConfigEqualRegions =
    util::Config("partitioner", "equal_regions") |
    util::Config("halo", 1);
  const auto meshConfigCubedSphere =
    util::Config("partitioner", "cubedsphere") |
    util::Config("halo", 1);

  // Set mesh generator.
  const auto meshGenEqualRegions = MeshGenerator("cubedsphere", meshConfigEqualRegions);
  const auto meshGenCubedSphere = MeshGenerator("cubedsphere", meshConfigCubedSphere);

  // Set mesh
  const auto meshEqualRegions = meshGenEqualRegions.generate(grid);
  const auto meshCubedSphere = meshGenCubedSphere.generate(grid);

  // Set functionspace.
  const auto equalRegionsCellColumns =
    functionspace::CubedSphereCellColumns(meshEqualRegions);
  const auto cubedSphereCellColumns =
    functionspace::CubedSphereCellColumns(meshCubedSphere);
  const auto equalRegionsNodeColumns =
    functionspace::CubedSphereNodeColumns(meshEqualRegions);
  const auto cubedSphereNodeColumns =
    functionspace::CubedSphereNodeColumns(meshCubedSphere);

  // test functionspaces.
  SECTION("CellColumns: equal_regions") {
    testFunctionSpace(equalRegionsCellColumns);
  }
  SECTION("CellColumns: cubedsphere") {
    testFunctionSpace(cubedSphereCellColumns);
  }
  SECTION("NodeColumns: equal_regions") {
    testFunctionSpace(equalRegionsNodeColumns);
  }
  SECTION("NodeColumns: cubedsphere") {
    testFunctionSpace(cubedSphereNodeColumns);
  }


}



}  // namespace test
}  // namespace atlas

int main( int argc, char** argv ) {
  return atlas::test::run( argc, argv );
}
