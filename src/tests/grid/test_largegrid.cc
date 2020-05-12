#include <algorithm>
#include <vector>

#include "atlas/array.h"
#include "atlas/functionspace.h"
#include "atlas/grid.h"
#include "atlas/library/Library.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/util/Config.h"
#include "atlas/util/UnitSphere.h"
#include "eckit/types/FloatCompare.h"
#include "tests/AtlasTestEnvironment.h"

using namespace atlas::util;
using namespace atlas::grid;
using namespace atlas;


namespace atlas {
namespace test {
CASE( "largeGrid" ) {



  atlas::StructuredGrid grid (atlas::util::Config ("nx", 40000) | atlas::util::Config ("ny", 20000) 
                            | atlas::util::Config ("type", "regular_lonlat") 
                            | atlas::util::Config ("domain", atlas::util::Config ("type", "global") 
                                                 | atlas::util::Config ("units", "degrees")
                                                 | atlas::util::Config ("west", -180.)));
 
                       
  auto lonlat = grid.lonlat (0, 0);

  printf (" %12.2f, %12.2f\n", lonlat.lon (), lonlat.lat ());


}
}  // namespace test
}  // namespace atlas


int main( int argc, char* argv[] ) {
    return atlas::test::run( argc, argv );
}
