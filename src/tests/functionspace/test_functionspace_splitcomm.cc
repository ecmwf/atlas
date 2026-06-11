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
#include "atlas/functionspace/PointCloud.h"
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
    return Grid(color() == 0 ? "O32" : "N32" );
}

std::string expected_checksum() {
    static std::string result = [&]() {
        if (grid().name()=="O32") {
            return "6408";
        }
        else if (grid().name()=="N32") {
            return "ca85";
        }
        else {
            return "unknown";
        }
    }();
    return result;
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
    auto g = array::make_view<gidx_t,1>(fs.global_index());

    if (functionspace::BlockStructuredColumns fsb{fs}) {
        auto f = array::make_view<double,2>(field);
        auto g = array::make_view<gidx_t,1>(fs.global_index());
        for( idx_t jblk=0; jblk<fsb.nblks(); ++jblk ) {
            auto block = fsb.block(jblk);
            for( idx_t jlane=0; jlane<block.size(); ++jlane ) {
                f(jblk,jlane) = g(block.index(jlane));
            }
        }
        return;
    }
    auto f = array::make_view<double,1>(field);
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

    Log::error() << "fs.part() = " << fs.part() << " fs.nb_parts() = " << fs.nb_parts() << " grid = " << fs.grid().name() << " checksum = " << checksum << std::endl;

    // Output
    output::Gmsh gmsh(mesh.grid().name()+".msh");
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

    Log::error() << "fs.part() = " << fs.part() << " fs.nb_parts() = " << fs.nb_parts() << " grid = " << fs.grid().name() << " checksum = " << checksum << std::endl;
}

CASE("test FunctionSpace BlockStructuredColumns") {
    Fixture fixture;

    auto fs   = functionspace::BlockStructuredColumns(grid(),atlas::option::halo(1)|option::mpi_split_comm());
    EXPECT_EQUAL(fs.part(),mpi::comm("split").rank());
    EXPECT_EQUAL(fs.nb_parts(),mpi::comm("split").size());

    auto field  = fs.createField<double>();
    field_init(field);

    // HaloExchange
    // field.haloExchange();
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

    Log::error() << "fs.part() = " << fs.part() << " fs.nb_parts() = " << fs.nb_parts() << " grid = " << fs.grid().name() << " checksum = " << checksum << std::endl;
}

//-----------------------------------------------------------------------------

CASE("test FunctionSpace PointCloud") {
    Fixture fixture;

    auto fs   = functionspace::PointCloud(grid(),util::Config("halo_radius",400*1000)|option::mpi_split_comm());
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
    // auto checksum = fs.checksum(field);
    // EXPECT_EQ(checksum, expected_checksum());
}

//-----------------------------------------------------------------------------

CASE("test FunctionSpace StructuredColumns with MatchingPartitioner") {
    Fixture fixture;
    auto g = grid();
    auto fs_A = functionspace::StructuredColumns(g, option::mpi_split_comm());
    auto fs_B = functionspace::StructuredColumns(g, grid::MatchingPartitioner(fs_A), option::mpi_split_comm());
    fs_A.polygon().outputPythonScript("fs_A_polygons.py");
    fs_B.polygon().outputPythonScript("fs_B_polygons.py");
}

//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
