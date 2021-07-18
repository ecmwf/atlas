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
      const auto grid = atlas::Grid("CS-EA-2");

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


    CASE("cubedsphere_tile_anyXYToFundamentalXY_test") {

      auto tileConfig1 = atlas::util::Config("type", "cubedsphere_lfric");
      auto lfricTiles = atlas::CubedSphereTiles(tileConfig1);

      // create a nodal cubed-sphere grid and check that no point are changed by
      // iterating through points.

      int resolution( 2 );
      std::vector<std::string> grid_names{"CS-LFR-" + std::to_string( resolution ),
                                         };
      Grid grid{grid_names[0]};

      int jn{0};
      for ( auto crd : grid.xy() ) {
          atlas::PointXY initialXY{crd[XX], crd[YY]};
          atlas::PointXY finalXY = lfricTiles.anyXYToFundamentalXY(initialXY);
          EXPECT_APPROX_EQ(initialXY, finalXY);
          ++jn;
      }

      std::vector<std::pair<double, double>> gettingXY{
          {180, -45}, {225.0,-45.0}
      };

      std::vector<std::pair<double, double>> expectedXY{
          {270.,-45.0}, {225.0, -45.0},  {180,-45.0},  {135.0,-45.0}, {90.0,-45.0},   {45.0, -45.0},   {270.0,-45.0}, {315.0,-45.0}, {270.0,-45.0},
          {315.,-45.0},  {45.0, -90.0},  {135,-45.0},  {45.0, -90.0}, {135.0, -45.0}, {45.0, -90.0}, {315.0,-45.0}, {45.0, -90.0}, {315.,-45.0},
          {0.0, -45.0},  {45.0, -45.0}, {90.0,-45.0},  {135.0,-45.0}, {180.0, -45.0}, {225.0,-45.0}, {270.0,-45.0}, {315.0,-45.0}, {0.0,-45.0},
          {0.0,   0.0},  {45.0,  0.0},  {90.0,   0.0}, {135.0,  0.0},   {180.0, 0.0},  {225.0, 0.0},  {270.0, 0.0}, {315.0,  0.0}, {0.0, 0.0},
          {0.0,   45.0}, {45.0,  45.0}, {90.0,  45.0},                                                                             {0.0, 45.0},
          {0.0,   90.0}, {45.0,  90.0}, {90.0,  90.0},                                                                             {0.0, 90.0},
          {0.0,  135.0}, {45.0, 135.0}, {90.0, 135.0},                                                                             {0.0,  135.0},
          {180.0,  0.0}, {225.0,  0.0}, {270.0,  0.0}, {315.0,  0.0}, {0.0,   0.0}, {45.0,  0.0},  {90.0,   0.0},  {135.0, 0.0},  {180.0, 0.0},
          {180.0, -45.0},{225.0,-45.0}, {270.0,-45.0}, {315.0,-45.0}, {0.0, -45.0}, {45.0, -45.0}, {90.0, -45.0},  {135.0, -45.0},{180.0, -45.0} };

      // iterate through 9 points in x from [0,360]
      // iterate through 9 points in y from [-135, 225]
      // and check with expected output. 81 checks!


      std::vector<atlas::PointXY> startingXY{
          {0.0, 45.0}, {-45.0, 90}, {0.0, -45.0}, {0.0, -45.0} };

      size_t indx{0};

      for (idx_t t = 0; t < 1; ++t) {
          // for each tile we check whether the xy points are consistent Up, Right, Down, Left
          // of the tile.

          if (t == 0) {
              std::vector<atlas::PointXY> expectedXYUP{
                  {0.0, 45.0}, {45.0, 45.0},  {90.0, 45.0},
                  {0.0, 90.0}, {45.0, 90.0},  {90.0, 90.0},
                  {0.0,135.0}, {45.0,135.0},  {90.0,135.0},
                  {270.0,-45.0},{225.0,-45.0}, {180.0,-45.0},
                  {315.0,-45.0},{45.0, -90.0}, {135.0,-45.0},
                  {0.0,-45.0},{45.0, -45.0},  {90.0,-45.0}
              };


              // Up
              for (idx_t yIndx = 0 ; yIndx < 7 ; ++yIndx) {
                  for (idx_t xIndx = 0 ; xIndx < 3 ; ++xIndx) {
                      atlas::PointXY offsetXY{xIndx*45.0, yIndx*45.0};
                      atlas::PointXY initialXY = startingXY[t] + offsetXY;
                      atlas::PointXY middleXY = lfricTiles.anyXYToFundamentalXY(initialXY);
                      atlas::PointXY finalXY = lfricTiles.anyXYToFundamentalXY(middleXY);

                      std::cout << "xIndx yIndx InitialXY MiddleXY FinalXY " << xIndx  << " " << yIndx << "   " << initialXY.x() << " " << initialXY.y()
                                << "   " << middleXY.x() << " " << middleXY.y()
                                << "   " << finalXY.x() << " " << finalXY.y() << std::endl;

                      EXPECT_APPROX_EQ(middleXY, finalXY);
                      EXPECT_APPROX_EQ(middleXY, expectedXYUP[indx]);
                      indx += 1;

                  }
              }

          }
      }



    }

 
  }  // namespace test
}  // namespace atlas

int main( int argc, char** argv ) {
  return atlas::test::run( argc, argv );
}
