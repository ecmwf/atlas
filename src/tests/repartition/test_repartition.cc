

#include "tests/AtlasTestEnvironment.h"

namespace atlas {
  namespace test {

    CASE ("hello world") {
     std::cout << "Hello, world!" << std::endl;
    }

  }
}


int main( int argc, char** argv ) {
    return atlas::test::run( argc, argv );
}
