#include <string>
#include <vector>

#include "atlas/grid.h"
#include "atlas/mesh/HybridElements.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/meshgenerator.h"
#include "atlas/output/Gmsh.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/util/Config.h"

#include "tests/AtlasTestEnvironment.h"

namespace atlas {
namespace test {

void write_mesh_with_empty_partition(bool gather) {
    EXPECT_EQ(mpi::size(), 2);

    Grid grid("F4");
    std::vector<int> partitions(grid.size(), 0);
    grid::Distribution distribution(2, grid.size(), partitions.data());
    Mesh mesh = StructuredMeshGenerator().generate(grid, distribution);

    if (mpi::rank() == 0) {
        EXPECT(mesh.cells().size() > 0);
    }
    else {
        EXPECT_EQ(mesh.cells().size(), 0);
        EXPECT_EQ(mesh.edges().size(), 0);
    }

    const std::string filename = gather ? "test_gmsh_empty_partition_gather.msh"
                                        : "test_gmsh_empty_partition.msh";
    output::Gmsh(filename, util::Config("gather", gather) | util::Config("info", true)).write(mesh);
    mpi::comm().barrier();
}

CASE("test_gmsh_empty_partition") {
    write_mesh_with_empty_partition(false);
}

CASE("test_gmsh_empty_partition_gather") {
    write_mesh_with_empty_partition(true);
}

}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}