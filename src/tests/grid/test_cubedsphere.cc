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

#include "tests/AtlasTestEnvironment.h"

namespace atlas {
  namespace test {


    CASE("cubedsphere_grid_mesh_field_test") {

      // THIS IS TEMPORARY!
      // I expect this will be replaced by some more aggressive tests.

      // Set grid.
      const auto grid = atlas::Grid("CS-EA-24");

      atlas::Log::info() << grid->type() << std::endl;
      atlas::Log::info() << grid.size() << std::endl;


      // Set mesh.
      auto meshGen = atlas::MeshGenerator("cubedsphere");
      auto mesh = meshGen.generate(grid);

      // Set functionspace
      auto functionSpace = atlas::functionspace::NodeColumns(mesh,
        atlas::util::Config("levels", 1) | atlas::util::Config("periodic_points", true));

      // Set field
      auto field = functionSpace.createField<idx_t>(atlas::option::name("indices"));
      auto fieldView = atlas::array::make_view<idx_t, 2>( field );

      for (idx_t i = 0; i < fieldView.shape()[0]; ++i) {
        fieldView(i, 0) = i;
      }

      // Set gmsh config.
      auto gmshConfigXy = atlas::util::Config("coordinates", "xy");
      auto gmshConfigXyz = atlas::util::Config("coordinates", "xyz");
      auto gmshConfigLonLat = atlas::util::Config("coordinates", "lonlat");

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

  }
}

int main( int argc, char** argv ) {
  return atlas::test::run( argc, argv );
}
