/*
 * (C) British Crown Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "atlas/field/FieldSet.h"
#include "atlas/functionspace/StructuredColumns.h"
#include "atlas/grid.h"
#include "atlas/grid/Partitioner.h"
#include "atlas/meshgenerator.h"
#include "atlas/output/Gmsh.h"
#include "atlas/repartition/Repartition.h"
#include "atlas/util/Config.h"

#include "tests/AtlasTestEnvironment.h"

namespace atlas {
  namespace test {


    // Define test pattern for grid.
    double testPattern(
      const double lambda, const double phi,
      const idx_t field, const idx_t level) {

      return (1 + field) * std::cos(lambda * (1 + level) * M_PI / 180.)
        * std::cos(phi * (1 + level) * M_PI / 180.);
    }

    // Test repartitioner. Return true if test passed.
    // Output fields if gmshOutput == true.
    bool testStructColsToStructCols(
      const atlas::Grid& grid, const idx_t nFields, const idx_t nLevels,
      const atlas::grid::Partitioner& sourcePartitioner,
      const atlas::grid::Partitioner& targetPartitioner,
      const bool gmshOutput = false, const std::string& fileId = "") {


      // Set up function spaces.
      auto funcSpaceConfig = atlas::util::Config{};
      funcSpaceConfig.set("halo", 1);
      funcSpaceConfig.set("periodic_points", true);
      funcSpaceConfig.set("levels", nLevels);

      const auto sourceFunctionSpace = atlas::functionspace::StructuredColumns(
        grid, sourcePartitioner, funcSpaceConfig);

      const auto targetFunctionSpace = atlas::functionspace::StructuredColumns(
        grid, targetPartitioner, funcSpaceConfig);

      // Generate some field sets.
      auto sourceFieldSet = atlas::FieldSet{};
      for (idx_t field = 0; field < nFields; ++field){
        sourceFieldSet.add(sourceFunctionSpace.createField<double>(
          atlas::option::name("source_field_" + std::to_string(field))));
      }

      auto targetFieldSet = atlas::FieldSet{};
      for (idx_t field = 0; field < nFields; ++field){
        targetFieldSet.add(targetFunctionSpace.createField<double>(
          atlas::option::name("target_field_" + std::to_string(field))));
      }

      // Write some data to source fields.
      for (idx_t field = 0; field < sourceFieldSet.size(); ++field) {

        auto fieldView =
          atlas::array::make_view<double, 2>(sourceFieldSet[field]);

        for (idx_t j = sourceFunctionSpace.j_begin();
          j < sourceFunctionSpace.j_end(); ++j) {

          for (idx_t i = sourceFunctionSpace.i_begin(j);
            i < sourceFunctionSpace.i_end(j); ++i) {

            for (idx_t level = 0;
              level < sourceFunctionSpace.levels(); ++level) {

              // get lon and lat.
              const auto xy = sourceFunctionSpace.compute_xy(i, j);
              const auto lonLat = grid.projection().lonlat(xy);
              const auto lon = lonLat.lon();
              const auto lat = lonLat.lat();
              const auto f = testPattern(lon, lat, field, level);

              // write f to field.
              const auto iNode = sourceFunctionSpace.index(i, j);
              fieldView(iNode, level) = f;

            }
          }
        }
      }

      // Set up repartitioner.
      auto repart =
        atlas::Repartition(sourceFunctionSpace, targetFunctionSpace);

      // Execute repartitioner.
      repart.execute(sourceFieldSet, targetFieldSet);

      // Read and data from target fields.
      bool testPassed = true;

      for (idx_t field = 0; field < targetFieldSet.size(); ++field) {

        auto fieldView =
          atlas::array::make_view<double, 2>(targetFieldSet[field]);

        for (idx_t j = targetFunctionSpace.j_begin();
          j < targetFunctionSpace.j_end(); ++j) {

          for (idx_t i = targetFunctionSpace.i_begin(j);
            i < targetFunctionSpace.i_end(j); ++i) {

            for (idx_t level = 0;
              level < targetFunctionSpace.levels(); ++level) {

              // get lon and lat.
              const auto xy = targetFunctionSpace.compute_xy(i, j);
              const auto lonLat = grid.projection().lonlat(xy);
              const auto lon = lonLat.lon();
              const auto lat = lonLat.lat();
              const auto g = testPattern(lon, lat, field, level);

              // read f from field.
              const auto iNode = targetFunctionSpace.index(i, j);
              auto f = fieldView(iNode, level);

              // check that f is *exactly* equal to g;
              testPassed = testPassed && (f == g);

            }
          }
        }
      }

      // Write mesh and fields to file.
      if (gmshOutput) {

          // Generate meshes.
          const auto meshGen = atlas::MeshGenerator("structured");
          const auto sourceMesh = meshGen.generate(grid, sourcePartitioner);
          const auto targetMesh = meshGen.generate(grid, targetPartitioner);

          // Set gmsh config.
          auto gmshConfig = atlas::util::Config{};
          gmshConfig.set("ghost", "true");

          // Write meshes.
          const auto sourceMeshFile = fileId + "_source_mesh.msh";
          atlas::output::Gmsh(sourceMeshFile, gmshConfig).write(sourceMesh);

          const auto targetMeshFile = fileId + "_target_mesh.msh";
          atlas::output::Gmsh(targetMeshFile, gmshConfig).write(targetMesh);

          // Write fields.
          for (idx_t field = 0; field < sourceFieldSet.size(); ++field) {

            const auto sourceFieldFile =
              fileId + sourceFieldSet[field].name() + ".msh";

            atlas::output::Gmsh(sourceFieldFile, gmshConfig)
              .write(sourceFieldSet[field], sourceFunctionSpace);

            const auto targetFieldFile =
              fileId + targetFieldSet[field].name() + ".msh";

            atlas::output::Gmsh(targetFieldFile, gmshConfig)
              .write(targetFieldSet[field], targetFunctionSpace);

          }
      }

      return testPassed;
    }

    CASE ("StructuredColumns to StructuredColumns") {


      SECTION("lonlat: checkerboard to equal_regions") {

        // Set grid.
        idx_t nFields = 5;
        idx_t nLevels = 10;

        auto grid = atlas::Grid("L48x37");

        // Set partitioners.
        auto sourcePartitioner = atlas::grid::Partitioner("checkerboard");
        auto targetPartitioner = atlas::grid::Partitioner("equal_regions");

        // Check repartitioner.
        EXPECT(testStructColsToStructCols(grid, nFields, nLevels,
          sourcePartitioner, targetPartitioner));

        return;
      }

      SECTION("lonlat: equal_regions to checkerboard") {

        idx_t nFields = 5;
        idx_t nLevels = 10;

        // Set grid.
        auto grid = atlas::Grid("L48x37");

        // Set partitioners.
        auto sourcePartitioner = atlas::grid::Partitioner("equal_regions");
        auto targetPartitioner = atlas::grid::Partitioner("checkerboard");

        // Check repartitioner.
        EXPECT(testStructColsToStructCols(grid, nFields, nLevels,
          sourcePartitioner, targetPartitioner));

        return;
      }

      SECTION("gaussian: equal_regions to equal_bands") {

        idx_t nFields = 5;
        idx_t nLevels = 10;

        // Set grid.
        auto grid = atlas::Grid("O16");

        // Set partitioners.
        auto sourcePartitioner = atlas::grid::Partitioner("equal_regions");
        auto targetPartitioner = atlas::grid::Partitioner("equal_bands");

        // Check repartitioner.
        EXPECT(testStructColsToStructCols(grid, nFields, nLevels,
          sourcePartitioner, targetPartitioner));

        return;
      }

      SECTION("gaussian: equal_bands to equal_regions") {

        idx_t nFields = 5;
        idx_t nLevels = 10;

        // Set grid.
        auto grid = atlas::Grid("O16");

        // Set partitioners.
        auto sourcePartitioner = atlas::grid::Partitioner("equal_bands");
        auto targetPartitioner = atlas::grid::Partitioner("equal_regions");

        // Check repartitioner.
        EXPECT(testStructColsToStructCols(grid, nFields, nLevels,
          sourcePartitioner, targetPartitioner));

        return;
      }

      SECTION("gaussian: gmsh output") {

        idx_t nFields = 1;
        idx_t nLevels = 5;

        // Set up gaussian grid.
        auto grid = atlas::Grid("O16");

        // Set partitioners.
        auto sourcePartitioner = atlas::grid::Partitioner("equal_bands");
        auto targetPartitioner = atlas::grid::Partitioner("equal_regions");

        // Check repartitioner.
        EXPECT(testStructColsToStructCols(grid, nFields, nLevels,
          sourcePartitioner, targetPartitioner, true, grid.name()));

        return;
      }

    }

  }
}


int main( int argc, char** argv ) {
    return atlas::test::run( argc, argv );
}
