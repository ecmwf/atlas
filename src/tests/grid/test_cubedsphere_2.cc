#include "atlas/grid/CubedSphereGrid2.h"
#include "tests/AtlasTestEnvironment.h"

namespace atlas {
namespace test {
namespace {

CASE("cubed sphere grid instantiation") {
    const auto grid = CubedSphereGrid2(12);
    grid.printCubedSphere2();

}
}  // namespace
}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) { return atlas::test::run(argc, argv); }
