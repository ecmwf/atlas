/*
 * (C) Crown Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "atlas/array/MakeView.h"
#include "atlas/field/FieldSet.h"
#include "atlas/functionspace/CellColumns.h"
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/functionspace/StructuredColumns.h"
#include "atlas/grid.h"
#include "atlas/grid/CubedSphereGrid.h"
#include "atlas/grid/Partitioner.h"
#include "atlas/grid/Tiles.h"
#include "atlas/grid/detail/partitioner/CubedSpherePartitioner.h"
#include "atlas/interpolation.h"
#include "atlas/mesh.h"
#include "atlas/meshgenerator.h"
#include "atlas/meshgenerator/detail/cubedsphere/CubedSphereUtility.h"
#include "atlas/option.h"
#include "atlas/output/Gmsh.h"
#include "atlas/util/CoordinateEnums.h"
#include "atlas/util/function/VortexRollup.h"
#include "tests/AtlasTestEnvironment.h"

#include "eckit/mpi/Operation.h"

namespace atlas {
namespace test {

CASE("cubedsphere_mesh_jacobian_test") {
    using namespace meshgenerator::detail::cubedsphere;

    // Set grid an N2 grid with a halo size of 1.
    const auto grid = atlas::Grid("CS-LFR-C-2");

    // Set Jacobian
    const auto jacobian = NeighbourJacobian(CubedSphereGrid(grid));

    // Set vectors of known good outputs.
    const auto xyLocalKgoVec = std::vector<PointXY>{
        {0, -90},   {45, -90},  {90, -90},  {-45, -45}, {135, -45}, {-45, 0},   {135, 0},   {-45, 45},  {135, 45},
        {0, 90},    {45, 90},   {90, 90},   {90, -90},  {135, -90}, {180, -90}, {45, -45},  {225, -45}, {45, 0},
        {225, 0},   {45, 45},   {225, 45},  {90, 90},   {135, 90},  {180, 90},  {315, -45}, {315, 0},   {315, 45},
        {270, -90}, {270, 90},  {225, -90}, {225, 90},  {180, -90}, {180, 90},  {135, -45}, {135, 0},   {135, 45},
        {405, -45}, {405, 0},   {405, 45},  {360, -90}, {360, 90},  {315, -90}, {315, 90},  {270, -90}, {270, 90},
        {225, -45}, {225, 0},   {225, 45},  {0, 0},     {45, 0},    {90, 0},    {-45, 45},  {135, 45},  {-45, 90},
        {135, 90},  {-45, 135}, {135, 135}, {0, 180},   {45, 180},  {90, 180},  {-45, -45}, {-45, -90}, {-45, -135},
        {0, 0},     {0, -180},  {45, 0},    {45, -180}, {90, 0},    {90, -180}, {135, -45}, {135, -90}, {135, -135}};

    const auto xyGlobalKgoVec = std::vector<PointXY>{
        {315, -45}, {45, -90}, {135, -45}, {315, -45}, {135, -45}, {315, 0},   {135, 0},   {0, 90},    {90, 90},
        {0, 90},    {45, 90},  {90, 90},   {45, -45},  {45, -90},  {225, -45}, {45, -45},  {225, -45}, {45, 0},
        {225, 0},   {45, 45},  {45, 135},  {45, 45},   {45, 90},   {45, 135},  {315, -45}, {315, 0},   {0, 90},
        {315, -45}, {0, 90},   {45, -90},  {45, 90},   {135, -45}, {90, 90},   {135, -45}, {135, 0},   {90, 90},
        {45, -45},  {45, 0},   {45, 45},   {45, -45},  {45, 45},   {45, -90},  {45, 90},   {225, -45}, {45, 135},
        {225, -45}, {225, 0},  {45, 135},  {0, 0},     {45, 0},    {90, 0},    {0, 0},     {90, 0},    {315, 0},
        {135, 0},   {270, 0},  {180, 0},   {270, 0},   {225, 0},   {180, 0},   {0, 0},     {315, 0},   {270, 0},
        {0, 0},     {270, 0},  {45, 0},    {225, 0},   {90, 0},    {180, 0},   {90, 0},    {135, 0},   {180, 0}};

    const auto ijGlobalKgoVec = std::vector<PointIJ>{
        {0, 1}, {1, 1}, {1, 0}, {0, 1}, {1, 0}, {1, 1}, {1, 1}, {0, 1}, {2, 1}, {0, 1}, {1, 1}, {2, 1},
        {1, 0}, {1, 1}, {0, 1}, {1, 0}, {0, 1}, {1, 1}, {1, 1}, {1, 0}, {1, 2}, {1, 0}, {1, 1}, {1, 2},
        {0, 1}, {1, 1}, {0, 1}, {0, 1}, {0, 1}, {1, 1}, {1, 1}, {1, 0}, {2, 1}, {1, 0}, {1, 1}, {2, 1},
        {1, 0}, {1, 1}, {1, 0}, {1, 0}, {1, 0}, {1, 1}, {1, 1}, {0, 1}, {1, 2}, {0, 1}, {1, 1}, {1, 2},
        {0, 1}, {1, 1}, {0, 1}, {0, 1}, {0, 1}, {1, 1}, {1, 1}, {1, 2}, {1, 2}, {1, 2}, {1, 1}, {1, 2},
        {0, 1}, {1, 1}, {1, 2}, {0, 1}, {1, 2}, {1, 1}, {1, 1}, {0, 1}, {1, 2}, {0, 1}, {1, 1}, {1, 2}};

    const auto tKgoVec = std::vector<idx_t>{3, 5, 1, 3, 1, 3, 1, 4, 4, 4, 4, 4, 0, 5, 2, 0, 2, 0, 2, 4, 4, 4, 4, 4,
                                            3, 3, 4, 3, 4, 5, 4, 1, 4, 1, 1, 4, 0, 0, 4, 0, 4, 5, 4, 2, 4, 2, 2, 4,
                                            0, 0, 1, 0, 1, 3, 1, 3, 2, 3, 2, 2, 0, 3, 3, 0, 3, 0, 2, 1, 2, 1, 1, 2};

    // Set kgo iterators.
    auto xyLocalKgoIt  = xyLocalKgoVec.cbegin();
    auto xyGlobalKgoIt = xyGlobalKgoVec.cbegin();
    auto ijGlobalKgoIt = ijGlobalKgoVec.cbegin();
    auto tKgoIt        = tKgoVec.cbegin();


    // Play around with some grids.
    for (idx_t t = 0; t < 6; ++t) {
        for (idx_t j = -1; j < 4; ++j) {
            for (idx_t i = -1; i < 4; ++i) {
                // Set ij object.
                const auto ij = PointIJ(i, j);

                // Only look at halo values.
                if (!jacobian.ijCross(ij)) {
                    continue;
                }
                if (jacobian.ijInterior(ij)) {
                    continue;
                }

                // Get known good outputs.
                const auto xyLocalKgo  = *xyLocalKgoIt++;
                const auto xyGlobalKgo = *xyGlobalKgoIt++;
                const auto ijGlobalKgo = *ijGlobalKgoIt++;
                const auto tKgo        = *tKgoIt++;

                // Test known good output.
                const auto xyLocal   = jacobian.xy(ij, t);
                const auto txyGlobal = jacobian.xyLocalToGlobal(xyLocal, t);
                const auto tijGlobal = jacobian.ijLocalToGlobal(ij, t);

                EXPECT(xyLocal == xyLocalKgo);
                EXPECT(txyGlobal.xy() == xyGlobalKgo);
                EXPECT(txyGlobal.t() == tKgo);
                EXPECT(tijGlobal.ij() == ijGlobalKgo);
                EXPECT(tijGlobal.t() == tKgo);

                // Check xy and ij transforms are consistent.
                EXPECT(jacobian.ij(txyGlobal.xy(), txyGlobal.t()) == jacobian.ijLocalToGlobal(ij, t).ij());

                EXPECT(jacobian.xy(tijGlobal.ij(), tijGlobal.t()) == jacobian.xyLocalToGlobal(xyLocal, t).xy());
            }
        }
    }
}

void testHaloExchange(const std::string& gridStr, const std::string& partitionerStr, idx_t halo, bool output = true) {
    // Set grid.
    const auto grid = Grid(gridStr);

    // Set mesh config.
    const auto meshConfig = util::Config("partitioner", partitionerStr) | util::Config("halo", halo);

    // Set mesh generator.
    const auto meshGen = MeshGenerator("cubedsphere", meshConfig);

    // Set mesh
    const auto mesh = meshGen.generate(grid);

    // Set function space.
    const auto nodeColumns = functionspace::NodeColumns(mesh);
    const auto cellColumns = functionspace::CellColumns(mesh);

    // ---------------------------------------------------------------------------
    // Test node columns halo exchange.
    // ---------------------------------------------------------------------------

    Log::info() << "Starting node columns test." << std::endl;

    // make a test field.
    auto testField1 = nodeColumns.createField<double>(util::Config("name", "test field (node columns)"));

    // Make some field views.
    auto testView1  = array::make_view<double, 1>(testField1);
    auto lonLatView = array::make_view<double, 2>(nodeColumns.lonlat());
    auto ghostView  = array::make_view<int, 1>(nodeColumns.ghost());

    // Set non-ghost values.
    idx_t testFuncCallCount = 0;
    for (idx_t i = 0; i < nodeColumns.size(); ++i) {
        // Stop once we've reached ghost points.
        if (ghostView(i)) {
            break;
        }

        testView1(i) = 1.0 + util::function::vortex_rollup(lonLatView(i, LON), lonLatView(i, LAT), 1.0);
        ++testFuncCallCount;
    }

    // Make sure that some of the field values were ghosts.
    if (halo > 0) {
        EXPECT(testFuncCallCount < nodeColumns.size());
    }

    nodeColumns.haloExchange(testField1);

    // Check all values after halo exchange.
    double maxError = 0;
    for (idx_t i = 0; i < nodeColumns.size(); ++i) {
        const double testVal = 1.0 + util::function::vortex_rollup(lonLatView(i, LON), lonLatView(i, LAT), 1.0);
        maxError             = std::max(maxError, std::abs(testView1(i) - testVal));
    }

    mpi::comm().allReduceInPlace(maxError, eckit::mpi::Operation::Code::MAX);
    Log::info() << "Test field max error (NodeColumns) : " << maxError << std::scientific << std::endl;
    EXPECT(maxError < 1e-12);

    Log::info() << "Passed node columns test." << std::endl;

    // ---------------------------------------------------------------------------
    // Test cell columns halo exchange.
    // ---------------------------------------------------------------------------

    Log::info() << "Starting cell columns test." << std::endl;

    // make a test field.
    auto testField2 = cellColumns.createField<double>(util::Config("name", "test field (cell columns)"));

    // Make some field views.
    auto testView2 = array::make_view<double, 1>(testField2);
    auto haloView  = array::make_view<int, 1>(cellColumns.mesh().cells().halo());
    lonLatView     = array::make_view<double, 2>(cellColumns.mesh().cells().field("lonlat"));

    // Set non-halo values.
    testFuncCallCount = 0;
    for (idx_t i = 0; i < cellColumns.size(); ++i) {
        if (haloView(i)) {
            break;
        }

        testView2(i) = 1.0 + util::function::vortex_rollup(lonLatView(i, LON), lonLatView(i, LAT), 1.0);
        ++testFuncCallCount;
    }

    // Make sure that some of the field values were ghosts.
    if (halo > 0) {
        EXPECT(testFuncCallCount < cellColumns.size());
    }

    cellColumns.haloExchange(testField2);

    // Check all values after halo exchange.
    maxError = 0;
    for (idx_t i = 0; i < cellColumns.size(); ++i) {
        // Test field and test function should be the same.
        const double testVal = 1.0 + util::function::vortex_rollup(lonLatView(i, LON), lonLatView(i, LAT), 1.0);
        maxError             = std::max(maxError, std::abs(testView2(i) - testVal));
    }

    mpi::comm().allReduceInPlace(maxError, eckit::mpi::Operation::Code::MAX);
    Log::info() << "Test field max error (CellColumns) : " << maxError << std::scientific << std::endl;
    EXPECT(maxError < 1e-12);

    Log::info() << "Passed cell columns test." << std::endl;


    if (output) {
        // Make a field set of useful diagnostic quantities.

        auto nodeFields = atlas::FieldSet{};
        nodeFields.add(mesh.nodes().xy());
        nodeFields.add(mesh.nodes().lonlat());
        nodeFields.add(mesh.nodes().ghost());
        nodeFields.add(mesh.nodes().halo());
        nodeFields.add(mesh.nodes().remote_index());
        nodeFields.add(mesh.nodes().partition());
        nodeFields.add(mesh.nodes().global_index());
        nodeFields.add(mesh.nodes().field("tij"));
        nodeFields.add(testField1);

        auto cellFields = FieldSet{};
        cellFields.add(mesh.cells().halo());
        cellFields.add(mesh.cells().remote_index());
        cellFields.add(mesh.cells().partition());
        cellFields.add(mesh.cells().global_index());
        cellFields.add(mesh.cells().field("tij"));
        cellFields.add(mesh.cells().field("xy"));
        cellFields.add(mesh.cells().field("lonlat"));
        cellFields.add(testField2);

        // Set gmsh config.
        const auto gmshConfigXy = util::Config("coordinates", "xy") | util::Config("ghost", true);

        const auto gmshConfigXyz = util::Config("coordinates", "xyz") | util::Config("ghost", true);

        // Set gmsh objects.
        const auto fileStr = gridStr + "_" + partitionerStr + "_halo" + std::to_string(halo);

        // Node columns output.
        auto gmshXy  = output::Gmsh(fileStr + "_NC_xy.msh", gmshConfigXy);
        auto gmshXyz = output::Gmsh(fileStr + "_NC_xyz.msh", gmshConfigXyz);

        gmshXy.write(mesh);
        gmshXy.write(nodeFields, nodeColumns);

        gmshXyz.write(mesh);
        gmshXyz.write(nodeFields, nodeColumns);

        // Cell columns output.
        gmshXy  = output::Gmsh(fileStr + "_CC_xy.msh", gmshConfigXy);
        gmshXyz = output::Gmsh(fileStr + "_CC_xyz.msh", gmshConfigXyz);

        gmshXy.write(mesh);
        gmshXy.write(cellFields, cellColumns);

        gmshXyz.write(mesh);
        gmshXyz.write(cellFields, cellColumns);
    }
}

CASE("cubedsphere_mesh_test") {
    SECTION("N12, halo = 0") {
        testHaloExchange("CS-LFR-C-12", "equal_regions", 0);
        testHaloExchange("CS-LFR-C-12", "cubedsphere", 0);
    }
    SECTION("N12, halo = 1") {
        testHaloExchange("CS-LFR-C-12", "equal_regions", 1);
        testHaloExchange("CS-LFR-C-12", "cubedsphere", 1);
    }
    SECTION("N12, halo = 3") {
        testHaloExchange("CS-LFR-C-12", "equal_regions", 3);
        testHaloExchange("CS-LFR-C-12", "cubedsphere", 3);
    }
    SECTION("Prime number mesh (N17)") {
        testHaloExchange("CS-LFR-C-17", "equal_regions", 1);
        testHaloExchange("CS-LFR-C-17", "cubedsphere", 1);
    }
}

CASE("cubedsphere_dual_mesh_test") {
    SECTION("NodeColumns") {
        // This test generates a test field on a structured gird, then interpolates
        // it to a cubed-sphere dual mesh field.
        // The following functionality is tested:
        //   * Generation of a cubed-sphere dual mesh.
        //   * Partitioning of cubed-sphere dual mesh with a matching mesh partitioner.
        //   * Interpolation from a StructuredColumns source to a cubed sphere dual mesh (NodeColumns) target.

        // Set number of levels for test field.
        const idx_t nLevels = 5;

        // Make a grid, mesh and functionspace for structured source grid.

        const auto sourceGrid          = Grid("O96");
        const auto sourceDistribution  = grid::Distribution(sourceGrid, grid::Partitioner("equal_bands"));
        const auto sourceMesh          = MeshGenerator("structured").generate(sourceGrid, sourceDistribution);
        const auto sourceFunctionSpace = functionspace::StructuredColumns(
            sourceGrid, sourceDistribution, util::Config("halo", 1) | util::Config("periodic_points", true));

        // Make and set a source field.
        auto sourceField            = sourceFunctionSpace.createField<double>(util::Config("name", "source field") |
                                                                   util::Config("levels", nLevels));
        auto sourceView             = array::make_view<double, 2>(sourceField);
        const auto sourceLonLatView = array::make_view<double, 2>(sourceFunctionSpace.lonlat());

        for (idx_t i = 0; i < sourceView.shape(0); ++i) {
            for (idx_t j = 0; j < sourceView.shape(1); ++j) {
                sourceView(i, j) =
                    1.0 + util::function::vortex_rollup(sourceLonLatView(i, LON), sourceLonLatView(i, LAT),
                                                        static_cast<double>(j) / (nLevels - 1));
            }
        }

        // Set matching mesh partitioner.
        const auto targetPartitioner = grid::MatchingPartitioner(sourceMesh);

        // Set target grid, mesh and functionspace.
        const auto targetGrid          = Grid("CS-LFR-C-48");
        const auto targetMesh          = MeshGenerator("cubedsphere_dual").generate(targetGrid, targetPartitioner);
        const auto targetFunctionSpace = functionspace::NodeColumns(targetMesh);
        auto targetField =
            targetFunctionSpace.createField<double>(util::Config("name", "targetField") | util::Config("levels", 5));

        // Perform interpolation.
        auto scheme = util::Config("type", "structured-linear2D") | util::Config("halo", 1);
        auto interp = Interpolation(scheme, sourceFunctionSpace, targetFunctionSpace);
        interp.execute(sourceField, targetField);
        targetField.haloExchange();

        // Test accuracy of interpolation.
        auto targetView             = array::make_view<double, 2>(targetField);
        const auto targetLonLatView = array::make_view<double, 2>(targetFunctionSpace.lonlat());
        double maxError             = 0.;
        for (idx_t i = 0; i < targetView.shape(0); ++i) {
            for (idx_t j = 0; j < targetView.shape(1); ++j) {
                double referenceVal =
                    1.0 + util::function::vortex_rollup(targetLonLatView(i, LON), targetLonLatView(i, LAT),
                                                        static_cast<double>(j) / (nLevels - 1));

                double relativeError = std::abs((targetView(i, j) - referenceVal) / referenceVal);
                maxError             = std::max(maxError, relativeError);
            }
        }
        mpi::comm().allReduceInPlace(maxError, eckit::mpi::Operation::Code::MAX);
        Log::info() << "Max interpolation error = " + std::to_string(maxError) << std::endl;

        // Maximum interpolation error should be less than 1%.
        EXPECT(maxError < 0.01);

        // gmsh output.
        const auto gmshConfigXy =
            util::Config("coordinates", "xy") | util::Config("ghost", true) | util::Config("info", true);
        const auto gmshConfigXyz =
            util::Config("coordinates", "xyz") | util::Config("ghost", true) | util::Config("info", true);
        auto gmshXy  = output::Gmsh("dual_nodes_xy.msh", gmshConfigXy);
        auto gmshXyz = output::Gmsh("dual_nodes_xyz.msh", gmshConfigXyz);
        gmshXy.write(targetMesh);
        gmshXyz.write(targetMesh);
        gmshXy.write(targetField, targetFunctionSpace);
        gmshXyz.write(targetField, targetFunctionSpace);
    }

    SECTION("CellColumns") {
        // Simple test to check halo exchange and output mesh.

        // Set number of levels for test field.
        const idx_t nLevels = 5;

        // Set grid, mesh and functionspace.
        const auto grid = Grid("CS-LFR-C-48");
        const auto mesh =
            MeshGenerator("cubedsphere_dual", util::Config("partitioner", "equal_regions")).generate(grid);
        const auto functionSpace = functionspace::CellColumns(mesh);
        auto field =
            functionSpace.createField<double>(util::Config("name", "targetField") | util::Config("levels", nLevels));

        // Make views.
        auto fieldView           = array::make_view<double, 2>(field);
        const auto lonLatView    = array::make_view<double, 2>(functionSpace.lonlat());
        const auto remoteIdxView = array::make_indexview<idx_t, 1>(functionSpace.cells().remote_index());
        const auto partView      = array::make_view<int, 1>(functionSpace.cells().partition());

        // Set test field on owned cells.
        for (idx_t i = 0; i < fieldView.shape(0); ++i) {
            // Break once we've got to ghost cells (cells mirrored from ghost nodes).
            if (partView(i) != mpi::rank() || remoteIdxView(i) != i) {
                break;
            }

            for (idx_t j = 0; j < fieldView.shape(1); ++j) {
                fieldView(i, j) = 1.0 + util::function::vortex_rollup(lonLatView(i, LON), lonLatView(i, LAT),
                                                                      static_cast<double>(j) / (nLevels - 1));
            }
        }

        // Perform halo exchange.
        field.haloExchange();

        // Check test field on *all* cells.
        for (idx_t i = 0; i < fieldView.shape(0); ++i) {
            for (idx_t j = 0; j < fieldView.shape(1); ++j) {
                EXPECT_APPROX_EQ(fieldView(i, j),
                                 1.0 + util::function::vortex_rollup(lonLatView(i, LON), lonLatView(i, LAT),
                                                                     static_cast<double>(j) / (nLevels - 1)));
            }
        }

        // gmsh output.
        const auto gmshConfigXy =
            util::Config("coordinates", "xy") | util::Config("ghost", true) | util::Config("info", true);
        const auto gmshConfigXyz =
            util::Config("coordinates", "xyz") | util::Config("ghost", true) | util::Config("info", true);
        auto gmshXy  = output::Gmsh("dual_cells_xy.msh", gmshConfigXy);
        auto gmshXyz = output::Gmsh("dual_cells_xyz.msh", gmshConfigXyz);
        gmshXy.write(mesh);
        gmshXyz.write(mesh);
        gmshXy.write(FieldSet{field}, functionSpace);
        gmshXyz.write(FieldSet{field}, functionSpace);
    }
}


}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
