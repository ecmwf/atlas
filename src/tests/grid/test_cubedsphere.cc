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


    CASE("cubedsphere_tile_constructor_test") {

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


    CASE("cubedsphere_tileCubePeriodicity_test") {

      auto tileConfig1 = atlas::util::Config("type", "cubedsphere_lfric");
      auto lfricTiles = atlas::CubedSphereTiles(tileConfig1);

      // create a nodal cubed-sphere grid and check that no point are changed by
      // iterating through points.

      int resolution( 2 );
      std::vector<std::string> grid_names{"CS-LFR-L-" + std::to_string( resolution ),
                                         };
      Grid grid{grid_names[0]};

      int jn{0};
      for ( auto crd : grid.xy() ) {
          atlas::PointXY initialXY{crd[XX], crd[YY]};
          double xy[2] = {initialXY.x(), initialXY.y()};
          atlas::idx_t t = lfricTiles.tileFromXY(xy);
          atlas::PointXY finalXY = lfricTiles.tileCubePeriodicity(initialXY, t);
          EXPECT_APPROX_EQ(initialXY, finalXY);
          ++jn;
      }

      std::vector<atlas::PointXY>
              startingXYTile0 { {0., 315.}, {90., 315.}, {0., 225.}, {90., 225.},
                           {0., 135.}, {90., 135.}, {45., 0.}, {45., 90.},
                           {0., -45.}, {90., -45.}, {0, -135.}, {90,-135.},
                           {-90., 45.}, {180., 45.}, {270., 45.}, {360., 45.},
                           {-90.,-45.}, {180.,-45.}, {270.,-45.}, {360.,-45.}
                         };

      std::vector<atlas::PointXY>
              expectedXYTile0 { {0., -45.}, {90., -45.}, {270., -45.}, {180., -45.},
                           {0., 135.}, {90., 135.}, {45., 0.}, {45., 90.},
                           {0., -45.}, {90., -45.}, {270., -45.}, {180., -45.},
                           {0., 135.}, {90., 135.}, {0., 135.}, {0., 45.},
                           {270., -45.}, {180.,-45.}, {270.,-45.}, {0.,-45.}
                         };

      std::vector<atlas::PointXY>
              startingXYTile1 { {90., 315.}, {180., 315.}, {90., 225.}, {180., 225.},
                           {90., 135.}, {180., 135.}, {135., 0.}, {135., 90.},
                           {90., -45.}, {180., -45.}, {90., -135.}, {180,-135.},
                           {0., 45.}, {270., 45.}, {360., 45.}, {450., 45.},
                           {0.,-45.}, {270.,-45.}, {360.,-45.}, {450.,-45.}
                         };

      std::vector<atlas::PointXY>
              expectedXYTile1 { {90., -45.}, {180., -45.}, {0., -45.}, {270., -45.},
                           {0., 45.}, {0., 135.}, {135., 0.}, {45., 90.},
                           {90., -45.}, {180., -45.}, {0., -45.}, {270., -45.},
                           {0.,  45.}, {0., 135.}, {0., 45.}, {90., 45.},
                           {0., -45.}, {270.,-45.}, {0.,-45.}, {90.,-45.}
                         };


      std::vector<atlas::PointXY>
              startingXYTile2 { {180., 315.}, {270., 315.}, {180., 225.}, {270., 225.},
                           {180., 135.}, {270., 135.}, {225., 0.}, {225., 90.},
                           {180., -45.}, {270., -45.}, {180., -135.}, {270,-135.},
                           {90., 45.}, {0., 45.}, {-90., 45.}, {180., 45.},
                           {90.,-45.}, {0.,-45.}, {90.,-45.}, {180.,-45.}
                         };

      std::vector<atlas::PointXY>
              expectedXYTile2{ {180., -45.}, {270., -45.}, {90., -45.}, {0., -45.},
                           {90., 45.}, {0., 45.}, {225., 0.}, {45., 90.},
                           {180., -45.}, {270., -45.}, {90., -45.}, {0., -45.},
                           {90.,  45.}, {0., 45.}, {0., 135.}, {90.0, 135.0},
                           {90., -45.}, {0.,-45.}, {90.,-45.}, {180.,-45.}
                         };

      std::vector<atlas::PointXY>
              startingXYTile3 { {270., 315.}, {360., 315.}, {270., 225.}, {360., 225.},
                           {270., 135.}, {360., 135.}, {315., 90.}, {315., -90.},
                           {270., -45.}, {360., -45.}, {270., -135.}, {360,-135.},
                           {180., 45.}, {90., 45.}, {360., 45.}, {270., 45.},
                           {180.,-45.}, {90.,-45.}, {360.,-45.}, {270.,-45.}
                         };

      std::vector<atlas::PointXY>
              expectedXYTile3{ {270., -45.}, {0., -45.}, {180., -45.}, {90., -45.},
                           {90., 135.}, {90., 45.}, {45., 90.}, {45., -90.},
                           {270., -45.}, {0., -45.}, {180., -45.}, {90., -45.},
                           {90.,  135.}, {90., 45.}, {0., 45.}, {0.0, 135.0},
                           {180., -45.}, {90.,-45.}, {0.,-45.}, {270.,-45.}
                         };


      std::vector<atlas::PointXY>
              startingXYTile4 { {0., 405.}, {90., 405.}, {0., 315.}, {90., 315.},
                           {0., 225.}, {90., 225.}, {0, -135.}, {45., -90.},
                           {0.,  45.}, {90.,  45.}, {0, -45.}, {90,-45.},
                           {-90., 135.}, {180., 135.}, {270., 135.}, {360., 135.},
                           {-90., 45.}, {180., 45.}, {270., 45.}, {360., 45.}
                         };

      std::vector<atlas::PointXY>
              expectedXYTile4 { {0., 45.}, {90., 45.}, {0., -45.}, {90., -45.},
                           {270., -45.}, {180., -45.}, {270., -45.}, {45., -90.},
                           {0., 45.}, {90., 45.}, {0., -45.}, {90., -45.},
                           {270., -45.}, {180., -45.}, {270., -45.}, {0., 135.},
                           {0., -45.}, {90., -45.}, {0., -45.}, {0., 45.}
                         };

      std::vector<atlas::PointXY>
              startingXYTile5 { {0.,225.}, {90., 225.}, {0., 135.}, {90., 135.},
                           {0., 45.}, {90., 45.}, {45., -90.}, {45., 0.},
                           {0., -135.}, {90., -135.}, {0, -225.}, {90,-225.},
                           {-90.,-45.}, {180.,-45.}, {270., -45.}, {360., -45.},
                           {-90.,-135.}, {180.,-135.}, {270.,-135.}, {360.,-135.}
                         };

      std::vector<atlas::PointXY>
              expectedXYTile5 { {270., -45.}, {180., -45.}, {0., 135.}, {90., 135.},
                           {0., 45.}, {90., 45.}, {45., -90.}, {45., 0.},
                           {270., -45.}, {180., -45.}, {0., 135.}, {90., 135.},
                           {0., 45.}, {90., 45.}, {0., 45.}, {0., -45.},
                           {0., 135.}, {90.,135.}, {0.,135.}, {270.,-45.}
                         };


      // testing tile 0
      for (atlas::idx_t t = 0; t < 6; t++) {
          std::vector<atlas::PointXY> startingXY;
          std::vector<atlas::PointXY> expectedXY;

          if (t == 0) {
              startingXY = startingXYTile0;
              expectedXY = expectedXYTile0;
          }
          if (t == 1) {
              startingXY = startingXYTile1;
              expectedXY = expectedXYTile1;
          }
          if (t == 2) {
              startingXY = startingXYTile2;
              expectedXY = expectedXYTile2;
          }
          if (t == 3) {
              startingXY = startingXYTile3;
              expectedXY = expectedXYTile3;
          }
          if (t == 4) {
              startingXY = startingXYTile4;
              expectedXY = expectedXYTile4;
          }
          if (t == 5) {
              startingXY = startingXYTile5;
              expectedXY = expectedXYTile5;
          }

          std::size_t jn{0};
          for (atlas::PointXY p : startingXY) {
              atlas::PointXY middleXY = lfricTiles.tileCubePeriodicity(p, t);
              atlas::PointXY finalXY = lfricTiles.tileCubePeriodicity(middleXY, 0);
              EXPECT_APPROX_EQ(middleXY, finalXY);
              EXPECT_APPROX_EQ(middleXY, expectedXY[jn]);
              ++jn;

          }

      }
    }

 
  }  // namespace test
}  // namespace atlas

int main( int argc, char** argv ) {
  return atlas::test::run( argc, argv );
}
