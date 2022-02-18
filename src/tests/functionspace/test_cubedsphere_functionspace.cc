/*
 * (C) Crown Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "atlas/array.h"
#include "atlas/field/FieldSet.h"
#include "atlas/functionspace/CellColumns.h"
#include "atlas/functionspace/CubedSphereColumns.h"
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/grid.h"
#include "atlas/grid/CubedSphereGrid.h"
#include "atlas/grid/Tiles.h"
#include "atlas/mesh.h"
#include "atlas/meshgenerator.h"
#include "atlas/meshgenerator/detail/cubedsphere/CubedSphereUtility.h"
#include "atlas/option.h"
#include "atlas/util/CoordinateEnums.h"
#include "atlas/util/function/VortexRollup.h"

#include "tests/AtlasTestEnvironment.h"

namespace atlas {
namespace test {

// Allow small differences in the last few bits of a double aprroximately equal to 1
constexpr double epsilon = std::numeric_limits<double>::epsilon() * 16;

template <typename BaseFunctionSpace>
void testFunctionSpace(const functionspace::CubedSphereColumns<BaseFunctionSpace>& functionspace) {
    // Make field.
    auto field     = functionspace.template createField<double>(util::Config("name", "test field"));
    auto fieldView = array::make_view<double, 1>(field);

    // Get view of lonlat.
    const auto lonLatView = array::make_view<double, 2>(functionspace.lonlat());

    // Get view of ghost field if NodeColumns, halo field if CellColumns.
    const auto ghostView = functionspace.type() == "NodeColumns"
                               ? array::make_view<int, 1>(functionspace.mesh().nodes().ghost())
                               : array::make_view<int, 1>(functionspace.mesh().nodes().halo());

    // Loop over all non halo elements of test field.
    idx_t testFuncCallCount = 0;
    functionspace.parallel_for([&](idx_t index, idx_t t, idx_t i, idx_t j) {
        // Make sure index matches ijt.
        EXPECT(index == functionspace.index(t, i, j));

        // Check that indices of "+" stencil are valid.
        const auto badIdx = functionspace.invalid_index();
        EXPECT(functionspace.index(t, i - 1, j) != badIdx);
        EXPECT(functionspace.index(t, i + 1, j) != badIdx);
        EXPECT(functionspace.index(t, i, j - 1) != badIdx);
        EXPECT(functionspace.index(t, i, j + 1) != badIdx);

        // Make sure we're avoiding halos.
        EXPECT(!ghostView(index));

        // Set field values.
        fieldView(index) = util::function::vortex_rollup(lonLatView(index, LON), lonLatView(index, LAT), 1.0);
        ++testFuncCallCount;
    });

    // Make sure call count is less than functionspace.size() as we skipped halos.
    EXPECT(testFuncCallCount < functionspace.size());

    // Perform halo exchange.
    functionspace.haloExchange(field);

    // Loop over elements including halo
    testFuncCallCount = 0;
    functionspace.parallel_for(util::Config("include_halo", true), [&](idx_t index, idx_t t, idx_t i, idx_t j) {
        // Make sure index matches ijt.
        EXPECT(index == functionspace.index(t, i, j));

        // Set field values.
        EXPECT_APPROX_EQ(fieldView(index),
                         util::function::vortex_rollup(lonLatView(index, LON), lonLatView(index, LAT), 1.0), epsilon);

        ++testFuncCallCount;
    });

    // Make sure call count is equal to functionspace.size().

    std::cout << testFuncCallCount << " " << functionspace.size() << std::endl;

    //EXPECT( testFuncCallCount == functionspace.size() );


    // Test SFINAE for parallel_for.
    // Suggestions for more inventive tests are welcome.

    idx_t nLevels   = 10;
    auto field2     = functionspace.template createField<double>(util::Config("name", "test field") |
                                                             util::Config("levels", nLevels));
    auto fieldView2 = array::make_view<double, 2>(field2);

    functionspace.parallel_for(util::Config("levels", nLevels), [&](idx_t index, idx_t t, idx_t i, idx_t j, idx_t k) {
        fieldView2(index, k) = t * i * j * k;
    });

    functionspace.parallel_for([&](idx_t index, idx_t t, idx_t i, idx_t j) {
        for (idx_t k = 0; k < nLevels; ++k) {
            EXPECT(static_cast<idx_t>(fieldView2(index, k)) == t * i * j * k);
        }
    });

    functionspace.parallel_for(util::Config("levels", nLevels),
                               [&](idx_t index, idx_t k) { fieldView2(index, k) = k; });

    functionspace.parallel_for([&](idx_t index) {
        for (idx_t k = 0; k < nLevels; ++k) {
            EXPECT(static_cast<idx_t>(fieldView2(index, k)) == k);
        }
    });
}

CASE("cubedsphere_mesh_functionspace") {
    // Set grid.
    const auto grid = Grid("CS-LFR-C-12");

    // Set mesh config.
    const auto meshConfigEqualRegions = util::Config("partitioner", "equal_regions") | util::Config("halo", 1);
    const auto meshConfigCubedSphere  = util::Config("partitioner", "cubedsphere") | util::Config("halo", 1);

    // Set mesh generator.
    const auto meshGenEqualRegions = MeshGenerator("cubedsphere", meshConfigEqualRegions);
    const auto meshGenCubedSphere  = MeshGenerator("cubedsphere", meshConfigCubedSphere);

    // Set dual mesh generator.
    const auto dualMeshGenEqualRegions =
        MeshGenerator("cubedsphere_dual", meshConfigEqualRegions);
    const auto dualMeshGenCubedSphere =
        MeshGenerator("cubedsphere_dual", meshConfigCubedSphere);

    // Set mesh
    const auto meshEqualRegions = meshGenEqualRegions.generate(grid);
    const auto meshCubedSphere  = meshGenCubedSphere.generate(grid);

    // Set dual mesh
    const auto dualMeshEqualRegions = dualMeshGenEqualRegions.generate(grid);
    const auto dualMeshCubedSphere  = dualMeshGenCubedSphere.generate(grid);

    // Set functionspace.
    const auto equalRegionsCellColumns     = functionspace::CubedSphereCellColumns(meshEqualRegions);
    const auto cubedSphereCellColumns      = functionspace::CubedSphereCellColumns(meshCubedSphere);
    const auto equalRegionsNodeColumns     = functionspace::CubedSphereNodeColumns(meshEqualRegions);
    const auto cubedSphereNodeColumns      = functionspace::CubedSphereNodeColumns(meshCubedSphere);
    const auto equalRegionsDualNodeColumns = functionspace::CubedSphereNodeColumns(dualMeshEqualRegions);
    const auto cubedSphereDualNodeColumns  = functionspace::CubedSphereNodeColumns(dualMeshCubedSphere);
    const auto equalRegionsDualCellColumns = functionspace::CubedSphereCellColumns(dualMeshEqualRegions);
    const auto cubedSphereDualCellColumns  = functionspace::CubedSphereCellColumns(dualMeshCubedSphere);


    // test functionspaces.
    SECTION("CellColumns: equal_regions") { testFunctionSpace(equalRegionsCellColumns); }
    SECTION("CellColumns: cubedsphere") { testFunctionSpace(cubedSphereCellColumns); }
    SECTION("NodeColumns: equal_regions") { testFunctionSpace(equalRegionsNodeColumns); }
    SECTION("NodeColumns: cubedsphere") { testFunctionSpace(cubedSphereNodeColumns); }
    SECTION("CellColumns: dual mesh, equal_regions") { testFunctionSpace(equalRegionsDualCellColumns); }
    SECTION("CellColumns: dual mesh, cubedsphere") { testFunctionSpace(cubedSphereDualCellColumns); }
    SECTION("NodeColumns: dual mesh, equal_regions") { testFunctionSpace(equalRegionsDualNodeColumns); }
    SECTION("NodeColumns: dual mesh, cubedsphere") { testFunctionSpace(cubedSphereDualNodeColumns); }
}

CASE("test copies and up/down casting") {
    // IMPORTANT: The internal structure should be cached and not be recomputed at all costs
    // This can be verified with ATLAS_DEBUG=1 and ATLAS_TRACE_REPORT=1

    // Set grid.
    const auto grid = Grid("CS-LFR-12");
    const auto mesh = MeshGenerator(grid.meshgenerator() | option::halo(1)).generate(grid);

    functionspace::CubedSphereCellColumns cs_cellcolumns;
    functionspace::CubedSphereNodeColumns cs_nodecolumns;

    ATLAS_TRACE_SCOPE("create functionspaces with cs structure") {
        cs_cellcolumns = functionspace::CubedSphereCellColumns(mesh);
        cs_nodecolumns = functionspace::CubedSphereNodeColumns(mesh);
    }

    ATLAS_TRACE_SCOPE("casting cellcolumns up/down") {
        const auto columns                               = functionspace::CellColumns(cs_cellcolumns);
        const auto functionspace                         = FunctionSpace(columns);
        const functionspace::CubedSphereCellColumns back = functionspace;
    }

    ATLAS_TRACE_SCOPE("casting cellcolumns up/down from mesh") {
        const auto columns                               = functionspace::CellColumns(mesh);
        const auto functionspace                         = FunctionSpace(columns);
        const functionspace::CubedSphereCellColumns back = functionspace;
    }

    ATLAS_TRACE_SCOPE("casting nodecolumns up/down") {
        const auto columns                               = functionspace::NodeColumns(cs_nodecolumns);
        const auto functionspace                         = FunctionSpace(columns);
        const functionspace::CubedSphereNodeColumns back = functionspace;
    }

    ATLAS_TRACE_SCOPE("casting nodecolumns up/down from mesh") {
        const auto columns                               = functionspace::NodeColumns(mesh);
        const auto functionspace                         = FunctionSpace(columns);
        const functionspace::CubedSphereNodeColumns back = functionspace;
    }
}

CASE("Cubed sphere primal-dual equivalence") {

    // Created a LFRic layout N12 cubedsphere grid.
    // This is a global set of grid points with lonlats located at cell centres.
    const auto grid = Grid("CS-LFR-12");

    // Create domain decomposed mesh.
    // We can generate two types of mesh, a primal and a dual.
    // The cell-centre (nodes) lonlats of the primal mesh are identical to the
    // node (cell-centre) lonlats of the dual mesh.
    const auto primalMesh = MeshGenerator("cubedsphere",
                                          util::Config("halo", 5) |
                                          util::Config("partitioner", "equal_regions")).generate(grid);
    const auto dualMesh = MeshGenerator("cubedsphere_dual",
                                        util::Config("halo", 4) |
                                        util::Config("partitioner", "equal_regions")).generate(grid);

    // Create cubed sphere function spaces (these have fancy features, such as
    // (t, i, j) indexing and parallel_for methods). The halo sizes of the primal
    // functionspaces are set to match that of the dual functionspaces.
    const auto primalNodes = functionspace::CubedSphereNodeColumns(primalMesh, util::Config("halo", 4));
    const auto primalCells = functionspace::CubedSphereCellColumns(primalMesh, util::Config("halo", 5));
    const auto dualNodes = functionspace::CubedSphereNodeColumns(dualMesh, util::Config("halo", 4));
    const auto dualCells = functionspace::CubedSphereCellColumns(dualMesh, util::Config("halo", 4));
    // Note, the functionspaces we are usually interested in are primalCells and
    // dualNodes. The others are there for completeness.

    // Field comparison function. Note that this upcasts everything to
    // to the base FunctionSpace class.
    const auto compareFields = [](const FunctionSpace& functionSpaceA,
                                  const FunctionSpace& functionSpaceB){

        // Check that function spaces are the same size.
        EXPECT_EQ(functionSpaceA.size(), functionSpaceB.size());

        // Make views to fields.
        const auto lonLatFieldA =
            array::make_view<double, 2>(functionSpaceA.lonlat());
        const auto lonLatFieldB =
            array::make_view<double, 2>(functionSpaceB.lonlat());
        const auto ghostFieldA =
            array::make_view<int, 1>(functionSpaceA.ghost());
        const auto ghostFieldB =
            array::make_view<int, 1>(functionSpaceB.ghost());

        // Loop over functionspaces.
        for (idx_t i = 0; i < functionSpaceA.size(); ++i) {
            EXPECT_EQ(lonLatFieldA(i, LON), lonLatFieldB(i, LON));
            EXPECT_EQ(lonLatFieldA(i, LAT), lonLatFieldB(i, LAT));
            EXPECT_EQ(ghostFieldA(i), ghostFieldB(i));
        }
    };

    // Check that primal cells and dual nodes are equivalent.
    compareFields(primalCells, dualNodes);

    // Check that dual cells and primal nodes are equivalent.
    compareFields(dualCells, primalNodes);

}

CASE("Variable halo size functionspaces (primal mesh)") {

    // Create a mesh with a large halo, and a few functionspaces with different
    // (smaller) halo sizes. These should create fields with a smaller memory
    // footprint.

    // Set grid.
    const auto grid = Grid("CS-LFR-C-12");

    // Set mesh config.
    const auto meshConfig = util::Config("partitioner", "equal_regions") | util::Config("halo", 3);

    // Set mesh.
    const auto mesh = MeshGenerator("cubedsphere", meshConfig).generate(grid);

    // Set functionspaces.
    const auto nodeColumns0 = functionspace::CubedSphereNodeColumns(mesh, util::Config("halo", 0));
    const auto nodeColumns1 = functionspace::CubedSphereNodeColumns(mesh, util::Config("halo", 1));
    const auto nodeColumns2 = functionspace::CubedSphereNodeColumns(mesh, util::Config("halo", 2));

    const auto cellColumns0 = functionspace::CubedSphereCellColumns(mesh, util::Config("halo", 0));
    const auto cellColumns1 = functionspace::CubedSphereCellColumns(mesh, util::Config("halo", 1));
    const auto cellColumns2 = functionspace::CubedSphereCellColumns(mesh, util::Config("halo", 2));

    // Check functionspace sizes.
    EXPECT(nodeColumns0.size() < nodeColumns1.size());
    EXPECT(nodeColumns1.size() < nodeColumns2.size());
    EXPECT(nodeColumns2.size() < mesh.nodes().size());
    EXPECT(cellColumns0.size() < cellColumns1.size());
    EXPECT(cellColumns1.size() < cellColumns2.size());
    EXPECT(cellColumns2.size() < mesh.cells().size());

    // Make sure size of owned cell data matches grid.
    auto checkSize = [&](idx_t sizeOwned){
        mpi::comm().allReduceInPlace(sizeOwned, eckit::mpi::Operation::SUM);
        EXPECT_EQ(sizeOwned, grid.size());
    };

    checkSize(cellColumns0.sizeOwned());
    checkSize(cellColumns1.sizeOwned());
    checkSize(cellColumns2.sizeOwned());

}

CASE("Variable halo size functionspaces (dual mesh)") {

    // Create a mesh with a large halo, and a few functionspaces with different
    // (smaller) halo sizes. These should create fields with a smaller memory
    // footprint.

    // Set grid.
    const auto grid = Grid("CS-LFR-C-12");

    // Set mesh config.
    const auto meshConfig = util::Config("partitioner", "equal_regions") | util::Config("halo", 3);

    // Set mesh.
    const auto mesh = MeshGenerator("cubedsphere_dual", meshConfig).generate(grid);

    // Set functionspaces.
    const auto nodeColumns0 = functionspace::CubedSphereNodeColumns(mesh, util::Config("halo", 0));
    const auto nodeColumns1 = functionspace::CubedSphereNodeColumns(mesh, util::Config("halo", 1));
    const auto nodeColumns2 = functionspace::CubedSphereNodeColumns(mesh, util::Config("halo", 2));

    const auto cellColumns0 = functionspace::CubedSphereCellColumns(mesh, util::Config("halo", 0));
    const auto cellColumns1 = functionspace::CubedSphereCellColumns(mesh, util::Config("halo", 1));
    const auto cellColumns2 = functionspace::CubedSphereCellColumns(mesh, util::Config("halo", 2));

    // Check functionspace sizes.
    EXPECT(nodeColumns0.size() < nodeColumns1.size());
    EXPECT(nodeColumns1.size() < nodeColumns2.size());
    EXPECT(nodeColumns2.size() < mesh.nodes().size());
    EXPECT(cellColumns0.size() < cellColumns1.size());
    EXPECT(cellColumns1.size() < cellColumns2.size());
    EXPECT(cellColumns2.size() < mesh.cells().size());

    // Make sure size of owned cell data matches grid.
    auto checkSize = [&](idx_t sizeOwned){
        mpi::comm().allReduceInPlace(sizeOwned, eckit::mpi::Operation::SUM);
        EXPECT_EQ(sizeOwned, grid.size());
    };

    checkSize(nodeColumns0.sizeOwned());
    checkSize(nodeColumns1.sizeOwned());
    checkSize(nodeColumns2.sizeOwned());

}

}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
