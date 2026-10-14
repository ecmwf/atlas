#include "atlas/functionspace/NodeColumns.h"
#include "atlas/grid/Distribution.h"
#include "atlas/grid/StructuredGrid.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/meshgenerator/MeshGenerator.h"
#include "atlas/grid/Partitioner.h"

#include "atlas/output/Gmsh.h"


#include "atlas/util/Config.h"
#include "tests/AtlasTestEnvironment.h"

using namespace atlas::mesh;

namespace atlas {
namespace test {

// ------------------------------------------------------------------

void test(std::string cs_gridname, std::string g_gridname) {
    const auto cs_grid        = Grid(cs_gridname);
    const auto cs_meshgen     = MeshGenerator{"cubedsphere_dual", util::Config("partitioner", "cubedsphere") | util::Config("halo", 1)};
    const auto cs_mesh        = cs_meshgen.generate(cs_grid);
    const auto cs_fs          = functionspace::NodeColumns(cs_mesh);
    const auto cs_partitioner = grid::MatchingMeshPartitioner(cs_fs.mesh(), option::type("cubedsphere"));

    const auto g_grid = StructuredGrid(g_gridname);
    const auto g_mesh = atlas::StructuredMeshGenerator().generate(g_grid, cs_partitioner);

    // For human verification
    output::Gmsh gmsh{g_gridname + "_with_" + cs_gridname + "_np" + std::to_string(cs_partitioner.nb_partitions()) + ".msh"};
    gmsh.write(g_mesh);
}

CASE("CS12 O12") {
    test("CS-LFR-12", "O12");
}

CASE("CS12 F12") {
    test("CS-LFR-12", "F12");
}

CASE("CS15 O15") {
    test("CS-LFR-15", "O15");
}

CASE("CS15 F15") {
    test("CS-LFR-15", "F15");
}

CASE("CS12 O8") {
    test("CS-LFR-12", "O8");
}

CASE("CS12 F8") {
    test("CS-LFR-12", "F8");
}

//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}