/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/array.h"
#include "atlas/grid.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/meshgenerator.h"
#include "atlas/option/Options.h"
#include "atlas/output/Gmsh.h"
#include "tests/AtlasTestEnvironment.h"

using namespace atlas::output;
using namespace atlas::meshgenerator;
using namespace atlas::grid;

namespace atlas {
namespace test {


//-----------------------------------------------------------------------------

namespace option {
    struct mpi_comm : public util::Config {
        mpi_comm(const atlas::mpi::Comm& comm) {
            set("mpi_comm",comm.name());
        }
        mpi_comm(std::string_view name) {
            set("mpi_comm",std::string(name));
        }
    };
}

CASE("test StructuredMeshGenerator directly") {
    auto& mpi_comm_world = mpi::comm("world");
    int color = mpi_comm_world.rank()%2;
    Grid grid;
    Mesh mesh;
    if (color == 0) {
        grid = Grid("O32");
    }
    else {
        grid = Grid("N32");
    }
    mpi_comm_world.split(color,"split");
    StructuredMeshGenerator meshgen{option::mpi_comm("split")};
    mesh = meshgen.generate(grid);
    EXPECT_EQUAL(mesh.nb_partitions(),2);
    EXPECT_EQUAL(mpi::comm().name(),"world");
    eckit::mpi::deleteComm("split");
}

CASE("test Mesh") {
    auto& mpi_comm_world = mpi::comm("world");
    int color = mpi_comm_world.rank()%2;
    Grid grid;
    if (color == 0) {
        grid = Grid("O32");
    }
    else {
        grid = Grid("N32");
    }
    mpi_comm_world.split(color,"split");
    auto mesh = Mesh(grid,option::mpi_comm("split"));
    EXPECT_EQUAL(mesh.nb_partitions(),2);
    EXPECT_EQUAL(mpi::comm().name(),"world");
    eckit::mpi::deleteComm("split");
    
}

//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
