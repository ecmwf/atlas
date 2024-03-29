/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/grid.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/meshgenerator.h"
#include "atlas/output/Gmsh.h"
#include "atlas/grid/Partitioner.h"

#include "tests/AtlasTestEnvironment.h"

namespace atlas {
namespace test {

//-----------------------------------------------------------------------------

namespace option {
    struct mpi_comm : public util::Config {
        mpi_comm(std::string_view name) {
            set("mpi_comm",std::string(name));
        }
    };
}

int color() {
    static int c = mpi::comm("world").rank()%2;
    return c;
}

Grid grid() {
    static Grid g (color() == 0 ? "O32" : "N32" );
    return g;
}

Grid grid_A() {
    return grid();
}

Grid grid_B() {
    static Grid g (color() == 0 ? "F64" : "O64" );
    return g;
}

Grid grid_CS() {
    static Grid g (color() == 0 ? "CS-LFR-C-8" : "CS-LFR-C-16" );
    return g;
}

Grid grid_Regular() {
    static Grid g (color() == 0 ? "S8" : "L16" );
    return g;
}

Grid grid_healpix() {
    static Grid g (color() == 0 ? "H8" : "H16" );
    return g;
}

struct Fixture {
    Fixture() {
        mpi::comm().split(color(),"split");
    }
    ~Fixture() {
        if (eckit::mpi::hasComm("split")) {
            eckit::mpi::deleteComm("split");
        }
    }
};

CASE("Partitioners") {
    Fixture fixture;
    SECTION("default") {
        auto partitioner = grid::Partitioner("equal_regions");
        EXPECT_EQ(partitioner.mpi_comm(),mpi::comm().name());
        auto distribution = partitioner.partition(Grid("O8"));
        EXPECT_EQ(distribution.nb_partitions(), mpi::comm().size());
    }
    SECTION("split") {
        auto partitioner = grid::Partitioner("equal_regions", option::mpi_comm("split"));
        EXPECT_EQ(partitioner.mpi_comm(),"split");
        auto distribution = partitioner.partition(Grid("O8"));
        EXPECT_EQ(distribution.nb_partitions(), mpi::comm("split").size());
    }
}

CASE("StructuredMeshGenerator") {
    Fixture fixture;

    MeshGenerator meshgen{"structured", option::mpi_comm("split")};
    Mesh mesh = meshgen.generate(grid());
    EXPECT_EQUAL(mesh.nb_parts(),mpi::comm("split").size());
    EXPECT_EQUAL(mesh.part(),mpi::comm("split").rank());
    EXPECT_EQUAL(mesh.mpi_comm(),"split");
    EXPECT_EQUAL(mpi::comm().name(),"world");
    output::Gmsh gmsh(grid().name()+"_1.msh");
    gmsh.write(mesh);

    // partitioning graph and polygon output
    EXPECT_NO_THROW(mesh.partitionGraph());
    EXPECT_NO_THROW(mesh.polygons());
    mesh.polygon().outputPythonScript(grid().name()+"_polygons_1.py");
}

CASE("RegularMeshGenerator") {
    Fixture fixture;

    MeshGenerator meshgen{"regular", option::mpi_comm("split")};
    Mesh mesh = meshgen.generate(grid_Regular());
    EXPECT_EQUAL(mesh.nb_parts(),mpi::comm("split").size());
    EXPECT_EQUAL(mesh.part(),mpi::comm("split").rank());
    EXPECT_EQUAL(mesh.mpi_comm(),"split");
    EXPECT_EQUAL(mpi::comm().name(),"world");
    output::Gmsh gmsh(grid().name()+"_regular.msh");
    gmsh.write(mesh);

    // partitioning graph and polygon output
    EXPECT_NO_THROW(mesh.partitionGraph());
    EXPECT_NO_THROW(mesh.polygons());
    mesh.polygon().outputPythonScript(grid().name()+"_regular_polygons_1.py");
}

CASE("HealpixMeshGenerator") {
    Fixture fixture;

    MeshGenerator meshgen{"healpix", option::mpi_comm("split")};
    Mesh mesh = meshgen.generate(grid_healpix());
    EXPECT_EQUAL(mesh.nb_parts(),mpi::comm("split").size());
    EXPECT_EQUAL(mesh.part(),mpi::comm("split").rank());
    EXPECT_EQUAL(mesh.mpi_comm(),"split");
    EXPECT_EQUAL(mpi::comm().name(),"world");
    output::Gmsh gmsh(grid().name()+"_regular.msh");
    gmsh.write(mesh);

    // partitioning graph and polygon output
    EXPECT_NO_THROW(mesh.partitionGraph());
    EXPECT_NO_THROW(mesh.polygons());
    mesh.polygon().outputPythonScript(grid().name()+"_regular_polygons_1.py");
}

CASE("DelaunayMeshGenerator") {
if( ATLAS_HAVE_TESSELATION ) {
    Fixture fixture;

    MeshGenerator meshgen{"delaunay", option::mpi_comm("split")};
    Mesh mesh = meshgen.generate(grid_CS());
    EXPECT_EQUAL(mesh.nb_parts(),mpi::comm("split").size());
    EXPECT_EQUAL(mesh.part(),mpi::comm("split").rank());
    EXPECT_EQUAL(mesh.mpi_comm(),"split");
    EXPECT_EQUAL(mpi::comm().name(),"world");
    output::Gmsh gmsh(grid_CS().name()+"_delaunay.msh",util::Config("coordinates","xyz"));
    gmsh.write(mesh);
}
else {
    Log::warning() << "Not testing DelaunayMeshGenerator as TESSELATION feature is OFF" << std::endl;
}
}

CASE("CubedSphereDualMeshGenerator") {
    Fixture fixture;

    MeshGenerator meshgen{"cubedsphere_dual", option::mpi_comm("split")};
    Mesh mesh = meshgen.generate(grid_CS());
    EXPECT_EQUAL(mesh.nb_parts(),mpi::comm("split").size());
    EXPECT_EQUAL(mesh.part(),mpi::comm("split").rank());
    EXPECT_EQUAL(mesh.mpi_comm(),"split");
    EXPECT_EQUAL(mpi::comm().name(),"world");
    output::Gmsh gmsh(grid_CS().name()+"_1.msh");
    gmsh.write(mesh);
}

CASE("CubedSphereMeshGenerator") {
    Fixture fixture;

    MeshGenerator meshgen{"cubedsphere", option::mpi_comm("split")};
    Mesh mesh = meshgen.generate(grid_CS());
    EXPECT_EQUAL(mesh.nb_parts(),mpi::comm("split").size());
    EXPECT_EQUAL(mesh.part(),mpi::comm("split").rank());
    EXPECT_EQUAL(mesh.mpi_comm(),"split");
    EXPECT_EQUAL(mpi::comm().name(),"world");
    output::Gmsh gmsh(grid_CS().name()+"_2.msh");
    gmsh.write(mesh);
}


CASE("Mesh constructor") {
    Fixture fixture;

    auto mesh = Mesh(grid(), option::mpi_comm("split"));
    EXPECT_EQUAL(mesh.nb_parts(),mpi::comm("split").size());
    EXPECT_EQUAL(mesh.part(),mpi::comm("split").rank());
    EXPECT_EQUAL(mesh.mpi_comm(),"split");
    EXPECT_EQUAL(mpi::comm().name(),"world");
    output::Gmsh gmsh(grid().name()+"_2.msh");
    gmsh.write(mesh);

    // partitioning graph and polygon output
    EXPECT_NO_THROW(mesh.partitionGraph());
    EXPECT_NO_THROW(mesh.polygons());
    mesh.polygon().outputPythonScript(grid().name()+"_polygons_2.py");
}

CASE("Mesh constructor with partitioner") {
    Fixture fixture;
    auto partitioner = grid::Partitioner("equal_regions", option::mpi_comm("split"));

    auto mesh = Mesh(grid(), partitioner);
    EXPECT_EQUAL(mesh.nb_parts(),mpi::comm("split").size());
    EXPECT_EQUAL(mesh.part(),mpi::comm("split").rank());
    EXPECT_EQUAL(mesh.mpi_comm(),"split");
    EXPECT_EQUAL(mpi::comm().name(),"world");
    output::Gmsh gmsh(grid().name()+"_2.msh");
    gmsh.write(mesh);

    // partitioning graph and polygon output
    EXPECT_NO_THROW(mesh.partitionGraph());
    EXPECT_NO_THROW(mesh.polygons());
    mesh.polygon().outputPythonScript(grid().name()+"_polygons_2.py");
}


CASE("MatchingPartitioner") {
    Fixture fixture;

    auto mesh_A = Mesh(grid_A(), option::mpi_comm("split"));
    auto mesh_B = Mesh(grid_B(), grid::MatchingPartitioner(mesh_A), option::mpi_comm("split"));

    output::Gmsh gmsh_B(grid_B().name()+"_3.msh");
    gmsh_B.write(mesh_B);

    // partitioning graph and polygon output
    EXPECT_NO_THROW(mesh_B.partitionGraph());
    EXPECT_NO_THROW(mesh_B.polygons());
    mesh_B.polygon().outputPythonScript(grid_B().name()+"_polygons_3.py");
}

}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
