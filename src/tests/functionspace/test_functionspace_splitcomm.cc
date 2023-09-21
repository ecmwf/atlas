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
#include "atlas/output/Gmsh.h"
#include "tests/AtlasTestEnvironment.h"
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/functionspace/StructuredColumns.h"
#include "atlas/functionspace/BlockStructuredColumns.h"
#include "atlas/grid/Partitioner.h"
#include "atlas/field/for_each.h"

namespace atlas {
namespace test {

//-----------------------------------------------------------------------------

namespace option {
    struct mpi_split_comm : public util::Config {
        mpi_split_comm() {
            set("mpi_comm","split");
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

std::string expected_checksum() {
    if (grid().name()=="O32") {
        return "75e913d400755a0d2782fc65e2035e97";
    }
    else if (grid().name()=="N32") {
        return "bcb344196d20becbb66f098d91f83abb";
    }
    else {
        return "unknown";
    }
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

void field_init(Field& field) {
    auto fs = field.functionspace();
    auto f = array::make_view<double,1>(field);
    auto g = array::make_view<gidx_t,1>(fs.global_index());
    for( idx_t j=0; j<f.size(); ++j ) {
        f(j) = g(j);
    }
}

CASE("test FunctionSpace NodeColumns") {
    Fixture fixture;

    auto mesh = Mesh(grid(), option::mpi_split_comm());
    auto fs = functionspace::NodeColumns(mesh,atlas::option::halo(1));
    EXPECT_EQUAL(fs.part(),mpi::comm("split").rank());
    EXPECT_EQUAL(fs.nb_parts(),mpi::comm("split").size());
    EXPECT_EQUAL(fs.mpi_comm(),"split");

    auto field  = fs.createField<double>();
    field_init(field);

    // HaloExchange
    field.haloExchange();
    // TODO CHECK

    // Gather
    auto fieldg = fs.createField<double>(atlas::option::global());
    fs.gather(field,fieldg);

    if (fieldg.size()) {
        idx_t g{0};
        field::for_each_value(fieldg,[&](double x) {
            EXPECT_EQ(++g,x);
        });
    }

    // Checksum
    auto checksum = fs.checksum(field);
    EXPECT_EQ(checksum, expected_checksum());

    // Output
    output::Gmsh gmsh(grid().name()+".msh");
    gmsh.write(mesh);
    gmsh.write(field);
}

CASE("test FunctionSpace StructuredColumns") {
    Fixture fixture;

    auto fs   = functionspace::StructuredColumns(grid(),atlas::option::halo(1)|option::mpi_split_comm());
    EXPECT_EQUAL(fs.part(),mpi::comm("split").rank());
    EXPECT_EQUAL(fs.nb_parts(),mpi::comm("split").size());

    auto field  = fs.createField<double>();
    field_init(field);

    // HaloExchange
    field.haloExchange();
    // TODO CHECK

    // Gather
    auto fieldg = fs.createField<double>(atlas::option::global());
    fs.gather(field,fieldg);

    if (fieldg.size()) {
        idx_t g{0};
        field::for_each_value(fieldg,[&](double x) {
            EXPECT_EQ(++g,x);
        });
    }

    // Checksum
    auto checksum = fs.checksum(field);
    EXPECT_EQ(checksum, expected_checksum());
}

CASE("test FunctionSpace BlockStructuredColumns") {
    Fixture fixture;

    auto fs   = functionspace::BlockStructuredColumns(grid(),atlas::option::halo(1)|option::mpi_split_comm());
    EXPECT_EQUAL(fs.part(),mpi::comm("split").rank());
    EXPECT_EQUAL(fs.nb_parts(),mpi::comm("split").size());

    auto field  = fs.createField<double>();
    // field_init(field);

    // HaloExchange
    // field.haloExchange();
    // TODO CHECK

    // Gather
    auto fieldg = fs.createField<double>(atlas::option::global());
    // fs.gather(field,fieldg);

    // if (fieldg.size()) {
    //     idx_t g{0};
    //     field::for_each_value(fieldg,[&](double x) {
    //         EXPECT_EQ(++g,x);
    //     });
    // }

    // // Checksum
    // auto checksum = fs.checksum(field);
    // EXPECT_EQ(checksum, expected_checksum());
}

//-----------------------------------------------------------------------------

CASE("test FunctionSpace StructuredColumns with MatchingPartitioner") {
    Fixture fixture;

    auto fs_A = functionspace::StructuredColumns(grid(), option::mpi_split_comm());
    auto fs_B = functionspace::StructuredColumns(grid(), grid::MatchingPartitioner(fs_A), option::mpi_split_comm());
    fs_A.polygon().outputPythonScript("fs_A_polygons.py");
    fs_B.polygon().outputPythonScript("fs_B_polygons.py");
}

//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
