/*
 * (C) Crown Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <cmath>

#include "atlas/field/FieldSet.h"
#include "atlas/functionspace/StructuredColumns.h"
#include "atlas/grid.h"
#include "atlas/grid/Partitioner.h"
#include "atlas/meshgenerator.h"
#include "atlas/output/Gmsh.h"
#include "atlas/redistribution/Redistribution.h"
#include "atlas/util/Config.h"

#include "tests/AtlasTestEnvironment.h"

namespace atlas {
namespace test {


int mpi_color() {
    static int c = mpi::comm("world").rank()%2;
    return c;
}

struct Fixture {
    Fixture() {
        mpi::comm().split(mpi_color(),"split");
    }
    ~Fixture() {
        if (eckit::mpi::hasComm("split")) {
            eckit::mpi::deleteComm("split");
        }
    }
};

// Define test pattern for grid.
template <typename T>
T testPattern(const double lambda, const double phi, const idx_t field, const idx_t level) {
    return static_cast<T>(100. * (1 + field) * std::cos(lambda * (1 + level) * M_PI / 180.) *
                          std::cos(phi * (1 + level) * M_PI / 180.));
}

// Define a default config for functionspaces.
atlas::util::Config funcSpaceDefaultConfig(const idx_t levels = 10, const idx_t halo = 1,
                                           const bool periodicPoints = true) {
    // Declare result.
    auto funcSpaceConfig = atlas::util::Config{};
    funcSpaceConfig.set("levels", levels);
    funcSpaceConfig.set("halo", halo);
    funcSpaceConfig.set("periodic_points", periodicPoints);

    return funcSpaceConfig;
}

// Test redistributer. Return true if test passed.
// Output fields if gmshOutput == true.
template <typename T>
bool testStructColsToStructCols(const atlas::Grid& grid, const idx_t nFields,
                                const atlas::grid::Partitioner& sourcePartitioner,
                                const atlas::grid::Partitioner& targetPartitioner,
                                const atlas::util::Config sourceFunctionSpaceConfig,
                                const atlas::util::Config targetFunctionSpaceConfig, const bool gmshOutput = false,
                                const std::string& fileId = "") {
    const auto sourceFunctionSpace =
        atlas::functionspace::StructuredColumns(grid, sourcePartitioner, sourceFunctionSpaceConfig);

    const auto targetFunctionSpace =
        atlas::functionspace::StructuredColumns(grid, targetPartitioner, targetFunctionSpaceConfig);

    // Generate some field sets.
    auto sourceFieldSet = atlas::FieldSet{};
    for (idx_t field = 0; field < nFields; ++field) {
        sourceFieldSet.add(
            sourceFunctionSpace.createField<T>(atlas::option::name("source_field_" + std::to_string(field))));
    }

    auto targetFieldSet = atlas::FieldSet{};
    for (idx_t field = 0; field < nFields; ++field) {
        targetFieldSet.add(
            targetFunctionSpace.createField<T>(atlas::option::name("target_field_" + std::to_string(field))));
    }

    // Write some data to source fields.
    for (idx_t field = 0; field < sourceFieldSet.size(); ++field) {
        auto fieldView = atlas::array::make_view<T, 2>(sourceFieldSet[field]);

        for (idx_t j = sourceFunctionSpace.j_begin(); j < sourceFunctionSpace.j_end(); ++j) {
            for (idx_t i = sourceFunctionSpace.i_begin(j); i < sourceFunctionSpace.i_end(j); ++i) {
                for (idx_t level = 0; level < sourceFunctionSpace.levels(); ++level) {
                    // get lon and lat.
                    const auto xy     = sourceFunctionSpace.compute_xy(i, j);
                    const auto lonLat = grid.projection().lonlat(xy);
                    const auto lon    = lonLat.lon();
                    const auto lat    = lonLat.lat();
                    const auto f      = testPattern<T>(lon, lat, field, level);

                    // write f to field.
                    const auto iNode        = sourceFunctionSpace.index(i, j);
                    fieldView(iNode, level) = f;
                }
            }
        }
    }

    // Set up redistributer.
    auto redist = atlas::Redistribution(sourceFunctionSpace, targetFunctionSpace,
                                        util::Config("type", "RedistributeStructuredColumns"));

    // Execute redistributer.
    redist.execute(sourceFieldSet, targetFieldSet);

    // Read and data from target fields.
    bool testPassed = true;

    for (idx_t field = 0; field < targetFieldSet.size(); ++field) {
        auto fieldView = atlas::array::make_view<T, 2>(targetFieldSet[field]);

        for (idx_t j = targetFunctionSpace.j_begin(); j < targetFunctionSpace.j_end(); ++j) {
            for (idx_t i = targetFunctionSpace.i_begin(j); i < targetFunctionSpace.i_end(j); ++i) {
                for (idx_t level = 0; level < targetFunctionSpace.levels(); ++level) {
                    // get lon and lat.
                    const auto xy     = targetFunctionSpace.compute_xy(i, j);
                    const auto lonLat = grid.projection().lonlat(xy);
                    const auto lon    = lonLat.lon();
                    const auto lat    = lonLat.lat();
                    const auto g      = testPattern<T>(lon, lat, field, level);

                    // read f from field.
                    const auto iNode = targetFunctionSpace.index(i, j);
                    auto f           = fieldView(iNode, level);

                    // check that f is *exactly* equal to g;
                    testPassed = testPassed && (f == g);
                }
            }
        }
    }

    // Write mesh and fields to file.
    if (gmshOutput) {
        // Generate meshes.
        const auto meshGen    = atlas::MeshGenerator("structured",util::Config("mpi_comm",sourceFunctionSpace.mpi_comm()));
        const auto sourceMesh = meshGen.generate(grid, sourcePartitioner);
        const auto targetMesh = meshGen.generate(grid, targetPartitioner);

        // Set gmsh config.
        auto gmshConfig = atlas::util::Config{};
        gmshConfig.set("ghost", "true");

        // Set source gmsh object.
        const auto sourceGmsh = atlas::output::Gmsh(fileId + "_source_mesh.msh", gmshConfig);


        // Set target gmsh object
        const auto targetGmsh = atlas::output::Gmsh(fileId + "_target_mesh.msh", gmshConfig);

        // Write gmsh
        sourceGmsh.write(sourceMesh);
        sourceGmsh.write(sourceFieldSet);
        targetGmsh.write(targetMesh);
        targetGmsh.write(targetFieldSet);
    }

    return testPassed;
}

CASE("Redistribute Structured Columns") {
    SECTION("lonlat: checkerboard to equal_regions") {
        // Set grid.
        idx_t nFields = 5;

        auto grid = atlas::Grid("L48x37");

        // Set partitioners.
        auto sourcePartitioner = atlas::grid::Partitioner("checkerboard");
        auto targetPartitioner = atlas::grid::Partitioner("equal_regions");

        // Check redistributer.
        EXPECT(testStructColsToStructCols<double>(grid, nFields, sourcePartitioner, targetPartitioner,
                                                  funcSpaceDefaultConfig(), funcSpaceDefaultConfig()));

        return;
    }

    SECTION("lonlat: equal_regions to checkerboard") {
        idx_t nFields = 5;

        // Set grid.
        auto grid = atlas::Grid("L48x37");

        // Set partitioners.
        auto sourcePartitioner = atlas::grid::Partitioner("equal_regions");
        auto targetPartitioner = atlas::grid::Partitioner("checkerboard");

        // Check redistributer.
        EXPECT(testStructColsToStructCols<double>(grid, nFields, sourcePartitioner, targetPartitioner,
                                                  funcSpaceDefaultConfig(), funcSpaceDefaultConfig()));

        return;
    }

    SECTION("gaussian: equal_regions to equal_bands") {
        idx_t nFields = 5;

        // Set grid.
        auto grid = atlas::Grid("O16");

        // Set partitioners.
        auto sourcePartitioner = atlas::grid::Partitioner("equal_regions");
        auto targetPartitioner = atlas::grid::Partitioner("equal_bands");

        EXPECT(testStructColsToStructCols<double>(grid, nFields, sourcePartitioner, targetPartitioner,
                                                  funcSpaceDefaultConfig(), funcSpaceDefaultConfig()));

        return;
    }

    SECTION("gaussian: equal_bands to equal_regions") {
        idx_t nFields = 5;

        // Set grid.
        auto grid = atlas::Grid("O16");

        // Set partitioners.
        auto sourcePartitioner = atlas::grid::Partitioner("equal_bands");
        auto targetPartitioner = atlas::grid::Partitioner("equal_regions");

        // Check redistributer.
        EXPECT(testStructColsToStructCols<double>(grid, nFields, sourcePartitioner, targetPartitioner,
                                                  funcSpaceDefaultConfig(), funcSpaceDefaultConfig()));

        return;
    }

    SECTION("gaussian: gmsh output") {
        idx_t nFields = 1;

        // Set up gaussian grid.
        auto grid = atlas::Grid("O16");

        // Set partitioners.
        auto sourcePartitioner = atlas::grid::Partitioner("equal_bands");
        auto targetPartitioner = atlas::grid::Partitioner("equal_regions");

        EXPECT(testStructColsToStructCols<double>(grid, nFields, sourcePartitioner, targetPartitioner,
                                                  funcSpaceDefaultConfig(), funcSpaceDefaultConfig(), true,
                                                  grid.name()));

        return;
    }

    SECTION("gaussian: <float> datatype") {
        idx_t nFields = 5;

        // Set grid.
        auto grid = atlas::Grid("O16");

        // Set partitioners.
        auto sourcePartitioner = atlas::grid::Partitioner("equal_bands");
        auto targetPartitioner = atlas::grid::Partitioner("equal_regions");

        EXPECT(testStructColsToStructCols<float>(grid, nFields, sourcePartitioner, targetPartitioner,
                                                 funcSpaceDefaultConfig(), funcSpaceDefaultConfig()));

        return;
    }

    SECTION("gaussian: <int> datatype") {
        idx_t nFields = 5;

        // Set grid.
        auto grid = atlas::Grid("O16");

        // Set partitioners.
        auto sourcePartitioner = atlas::grid::Partitioner("equal_bands");
        auto targetPartitioner = atlas::grid::Partitioner("equal_regions");

        EXPECT(testStructColsToStructCols<int>(grid, nFields, sourcePartitioner, targetPartitioner,
                                               funcSpaceDefaultConfig(), funcSpaceDefaultConfig()));

        return;
    }

    SECTION("gaussian: <long> datatype") {
        idx_t nFields = 5;

        // Set grid.
        auto grid = atlas::Grid("O16");

        // Set partitioners.
        auto sourcePartitioner = atlas::grid::Partitioner("equal_bands");
        auto targetPartitioner = atlas::grid::Partitioner("equal_regions");

        EXPECT(testStructColsToStructCols<long>(grid, nFields, sourcePartitioner, targetPartitioner,
                                                funcSpaceDefaultConfig(), funcSpaceDefaultConfig()));

        return;
    }

    SECTION("gaussian: mixed halo size I") {
        idx_t nFields = 5;

        // Set grid.
        auto grid = atlas::Grid("O16");

        // Set partitioners.
        auto sourcePartitioner = atlas::grid::Partitioner("equal_bands");
        auto targetPartitioner = atlas::grid::Partitioner("equal_regions");

        // set function space configs
        auto sourceFunctionSpaceConfig = funcSpaceDefaultConfig();
        sourceFunctionSpaceConfig.set("halo", 2);
        auto targetFunctionSpaceConfig = funcSpaceDefaultConfig();
        targetFunctionSpaceConfig.set("halo", 3);

        EXPECT(testStructColsToStructCols<double>(grid, nFields, sourcePartitioner, targetPartitioner,
                                                  sourceFunctionSpaceConfig, targetFunctionSpaceConfig));

        return;
    }

    SECTION("gaussian: mixed halo size II") {
        idx_t nFields = 5;

        // Set grid.
        auto grid = atlas::Grid("O16");

        // Set partitioners.
        auto sourcePartitioner = atlas::grid::Partitioner("equal_bands");
        auto targetPartitioner = atlas::grid::Partitioner("equal_regions");

        // set function space configs
        auto sourceFunctionSpaceConfig = funcSpaceDefaultConfig();
        sourceFunctionSpaceConfig.set("halo", 3);
        auto targetFunctionSpaceConfig = funcSpaceDefaultConfig();
        targetFunctionSpaceConfig.set("halo", 2);

        EXPECT(testStructColsToStructCols<double>(grid, nFields, sourcePartitioner, targetPartitioner,
                                                  sourceFunctionSpaceConfig, targetFunctionSpaceConfig));

        return;
    }
}

CASE("Redistribute Structured Columns with split comms") {
    Fixture fixture;
    SECTION("lonlat: checkerboard to equal_regions") {
        util::Config mpi_comm("mpi_comm","split");
        std::string id = std::to_string(mpi_color());

        // Set grid.
        idx_t nFields = 5;

        auto grid = atlas::Grid("L48x37");

        // Set partitioners.
        auto sourcePartitioner = atlas::grid::Partitioner(option::type("checkerboard")|mpi_comm);
        auto targetPartitioner = atlas::grid::Partitioner(option::type("equal_regions")|mpi_comm);

        // Check redistributer.
        EXPECT(testStructColsToStructCols<double>(grid, nFields, sourcePartitioner, targetPartitioner,
                                                  funcSpaceDefaultConfig(), funcSpaceDefaultConfig(),
                                                  true, id));
    }
}

}  // namespace test
}  // namespace atlas


int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
