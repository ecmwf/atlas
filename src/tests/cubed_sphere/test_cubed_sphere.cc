#include "atlas/functionspace/NodeColumns.h"
#include "atlas/grid.h"
#include "atlas/grid/Partitioner.h"
#include "atlas/mesh.h"
#include "atlas/meshgenerator.h"
#include "atlas/output/Gmsh.h"

#include "tests/AtlasTestEnvironment.h"

namespace atlas {
  namespace test {


    CASE("Cubed-sphere test.") {

      // Set grid.
      const auto grid = atlas::Grid("CS-EA-24");

      std::cout << grid->type() << std::endl;
      std::cout << grid.size() << std::endl;


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
