#include "atlas/functionspace/StructuredColumns.h"
#include "atlas/grid/Grid.h"
#include "atlas/grid/Partitioner.h"
#include "atlas/repartition/Repartition.h"
#include "atlas/util/Config.h"

#include "tests/AtlasTestEnvironment.h"

namespace atlas {
  namespace test {

    CASE ("hello world") {
      std::cout << "Hello, world!" << std::endl;

      // Build grids.
      atlas::util::Config gridOneConfig;
      gridOneConfig.set("type", "regular_lonlat");
      gridOneConfig.set("nx", 48);
      gridOneConfig.set("ny", 37);

      atlas::Grid gridOne(gridOneConfig);

      std::cout << "size " << gridOne.size() << std::endl;
      std::cout << "name " << gridOne.name() << std::endl;

      // Function spaces.
      atlas::util::Config funcSpaceOneConfig;
      funcSpaceOneConfig.set("halo", 1);
      funcSpaceOneConfig.set("periodic_points", true);
      funcSpaceOneConfig.set("levels", 1);

      atlas::functionspace::StructuredColumns
              funcSpaceOne(gridOne, atlas::grid::Partitioner("equal_regions"), funcSpaceOneConfig);

      atlas::functionspace::StructuredColumns
             funcSpaceTwo(gridOne, atlas::grid::Partitioner("checkerboard"), funcSpaceOneConfig);

      std::cout << "levels" << funcSpaceOne.levels() << std::endl;

      auto fieldOne = funcSpaceOne.createField<double>(atlas::option::name("field one"));
      auto fieldTwo = funcSpaceTwo.createField<double>(atlas::option::name("field two"));

      auto Repart = atlas::Repartition(funcSpaceOne, funcSpaceTwo);

      Repart.execute(fieldOne, fieldTwo);




    }

  }
}


int main( int argc, char** argv ) {
    return atlas::test::run( argc, argv );
}
