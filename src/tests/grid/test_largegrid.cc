#include <algorithm>
#include <vector>

#include <time.h>
#include <sys/resource.h>

#include "atlas/grid.h"
#include "atlas/library/Library.h"
#include "atlas/util/Config.h"
#include "tests/AtlasTestEnvironment.h"

namespace atlas {
namespace test {
CASE( "largeGrid" ) {

  time_t t0 = time (NULL);

  atlas::StructuredGrid grid (atlas::util::Config ("nx", 40000) | atlas::util::Config ("ny", 20000) 
                            | atlas::util::Config ("type", "regular_lonlat") 
                            | atlas::util::Config ("domain", atlas::util::Config ("type", "global") 
                                                 | atlas::util::Config ("units", "degrees")
                                                 | atlas::util::Config ("west", -180.)));
 
  EXPECT ((time (NULL)) - t0 < 10);

  struct rusage usage; 

  EXPECT (getrusage (RUSAGE_SELF, &usage) == 0);
  EXPECT (usage.ru_maxrss < 30000); // 30Mb

}
}  // namespace test
}  // namespace atlas


int main( int argc, char* argv[] ) {
    return atlas::test::run( argc, argv );
}
