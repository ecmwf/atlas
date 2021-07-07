/*
 * (C) Crown Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "atlas/array/MakeView.h"
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/grid.h"
#include "atlas/mesh.h"
#include "atlas/meshgenerator.h"
#include "atlas/output/Gmsh.h"
#include "atlas/grid/Tiles.h"
#include "atlas/option.h"

#include "tests/AtlasTestEnvironment.h"

namespace atlas {
  namespace test {

    CASE("cubedsphere_tile_test") {

      auto tileConfig1 = atlas::util::Config("type", "cubedsphere_lfric");
      auto lfricTiles = atlas::CubedSphereTiles(tileConfig1);
      EXPECT(lfricTiles.type() == "cubedsphere_lfric");

      auto tileConfig2 = atlas::util::Config("type", "cubedsphere_fv3");
      auto fv3Tiles = atlas::CubedSphereTiles(tileConfig2);
      EXPECT(fv3Tiles.type() == "cubedsphere_fv3");

      auto lfricTiles2 = atlas::CubedSphereTiles("cubedsphere_lfric");
      EXPECT(lfricTiles2.type() == "cubedsphere_lfric");

      auto fv3Tiles2 = atlas::CubedSphereTiles("cubedsphere_fv3");
      EXPECT(fv3Tiles.type() == "cubedsphere_fv3");

    }

    CASE("cubedsphere_grid_mesh_field_test") {

      // THIS IS TEMPORARY!
      // I expect this will be replaced by some more aggressive tests.

      // Set grid.
      const auto grid = atlas::Grid("CS-EA-L-2");

      atlas::Log::info() << grid->type() << std::endl;
      atlas::Log::info() << grid.size() << std::endl;


      // Set mesh.
      auto meshGen = atlas::MeshGenerator("cubedsphere");
      auto mesh = meshGen.generate(grid);

      // Set functionspace
      auto functionSpace = atlas::functionspace::NodeColumns(mesh);

      auto ghostIdx = mesh.nodes().metadata().get<std::vector<idx_t>>("ghost-global-idx");
      auto ownedIdx = mesh.nodes().metadata().get<std::vector<idx_t>>("owned-global-idx");

      // Print out ghost global indices with corresponding owned global indices
      auto ownedIdxIt = ownedIdx.begin();
      for (auto iGhost : ghostIdx) std::cout << iGhost << " " << *ownedIdxIt++ << std::endl;

      // Set field
      auto field = functionSpace.ghost();


      // Set gmsh config.
      auto gmshConfigXy = atlas::util::Config("coordinates", "xy") | atlas::util::Config("ghost", false);
      auto gmshConfigXyz = atlas::util::Config("coordinates", "xyz") | atlas::util::Config("ghost", false);
      auto gmshConfigLonLat = atlas::util::Config("coordinates", "lonlat") | atlas::util::Config("ghost", false);

      // Set source gmsh object.
      const auto gmshXy =
        atlas::output::Gmsh("cs_xy_mesh.msh", gmshConfigXy);
      const auto gmshXyz =
        atlas::output::Gmsh("cs_xyz_mesh.msh", gmshConfigXyz);
      const auto gmshLonLat =
        atlas::output::Gmsh("cs_lonlat_mesh.msh", gmshConfigLonLat);

      // Write gmsh.
      gmshXy.write(mesh);
      gmshXy.write(field);
      gmshXyz.write(mesh);
      gmshXyz.write(field);
      gmshLonLat.write(mesh);
      gmshLonLat.write(field);


    }

 
  }  // namespace test
}  // namespace atlas

int main( int argc, char** argv ) {
  return atlas::test::run( argc, argv );
}
