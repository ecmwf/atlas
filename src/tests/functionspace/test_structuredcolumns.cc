/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "eckit/log/Bytes.h"
#include "eckit/types/Types.h"

#include "atlas/array/ArrayView.h"
#include "atlas/array/MakeView.h"
#include "atlas/field/FieldSet.h"
#include "atlas/functionspace/BlockStructuredColumns.h"
#include "atlas/field/Field.h"
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/functionspace/StructuredColumns.h"
#include "atlas/grid/Partitioner.h"
#include "atlas/grid/StructuredGrid.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/meshgenerator.h"
#include "atlas/output/Gmsh.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/util/CoordinateEnums.h"
#include "atlas/util/MicroDeg.h"

#include "tests/AtlasTestEnvironment.h"

using namespace eckit;
using namespace atlas::functionspace;
using namespace atlas::util;

namespace atlas {
namespace test {

namespace {

template <typename Value>
void fill_structured_field(Field& field, Value first_value = Value{1}) {
    Value next_value = first_value;
    if (field.variables() && field.levels()) {
        auto value = array::make_view<Value, 3>(field);
        for (idx_t point = 0; point < field.shape(0); ++point) {
            for (idx_t jlev = 0; jlev < field.shape(1); ++jlev) {
                for (idx_t jvar = 0; jvar < field.shape(2); ++jvar) {
                    value(point, jlev, jvar) = next_value;
                    next_value += Value{1};
                }
            }
        }
    }
    else if (field.variables() || field.levels()) {
        auto value = array::make_view<Value, 2>(field);
        for (idx_t point = 0; point < field.shape(0); ++point) {
            for (idx_t entry = 0; entry < field.shape(1); ++entry) {
                value(point, entry) = next_value;
                next_value += Value{1};
            }
        }
    }
    else {
        auto value = array::make_view<Value, 1>(field);
        for (idx_t point = 0; point < field.shape(0); ++point) {
            value(point) = next_value;
            next_value += Value{1};
        }
    }
}

}  // namespace

//-----------------------------------------------------------------------------


CASE("ATLAS-295 (Github PR 32): Fix StructuredColumns when domain is shifted by 180 degrees (Github PR 32)") {
    auto grid1 = Grid("S80x40", GlobalDomain{-180});

    functionspace::StructuredColumns fs1(grid1, grid::Partitioner("checkerboard"),
                                         Config("halo", 1) | Config("periodic_points", true));

    auto grid2 = Grid{"S80x40"};

    functionspace::StructuredColumns fs2(grid2, grid::Partitioner("checkerboard"),
                                         Config("halo", 1) | Config("periodic_points", true));

    EXPECT_EQ(grid1.size(), grid2.size());
    EXPECT_EQ(fs1.size(), fs2.size());
    EXPECT_EQ(fs1.sizeOwned(), fs2.sizeOwned());
}


CASE("test_functionspace_StructuredColumns_no_halo") {
    size_t root          = 0;
    std::string gridname = eckit::Resource<std::string>("--grid", "O8");
    Grid grid(gridname);
    util::Config config;
    config.set("halo", 0);
    config.set("periodic_points", true);
    functionspace::StructuredColumns fs(grid, grid::Partitioner("equal_regions"), config);
    ATLAS_DEBUG_VAR(fs.size());
    ATLAS_DEBUG_VAR(eckit::Bytes(fs.footprint()));

    Field field     = fs.createField<double>(option::name("field"));
    Field field_glb = fs.createField<double>(option::name("field_global") | option::global(root));

    auto value     = array::make_view<double, 1>(field);
    auto value_glb = array::make_view<double, 1>(field_glb);

    value.assign(mpi::comm().rank());

    fs.gather(field, field_glb);

    Log::info() << "field checksum = " << fs.checksum(field) << std::endl;

    //  for( size_t j=0; j<value_glb.size(); ++j )
    //    Log::info() << value_glb(j) << " ";
    //  Log::info() << std::endl;

    if (mpi::comm().rank() == root && mpi::comm().size() == 5) {
        std::vector<double> check{
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3,
            3, 3, 3, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3,
            3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2,
            2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
            1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
            1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3,
            3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2,
            2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
            1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 1, 1, 1, 1,
            1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4,
            4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
            4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
            4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4};

        EXPECT(value_glb.size() == idx_t(check.size()));
        for (idx_t j = 0; j < value_glb.size(); ++j) {
            EXPECT(value_glb(j) == check[j]);
        }
    }
    ATLAS_TRACE_SCOPE("output gmsh") {
        output::Gmsh gmsh("structured.msh");

        gmsh.write(MeshGenerator("structured").generate(grid));
        gmsh.write(field);
    }
}

CASE("test_functionspace_StructuredColumns_halo with output") {
    ATLAS_DEBUG_VAR(mpi::comm().size());
    //  grid::StructuredGrid grid(
    //      grid::StructuredGrid::XSpace( {0.,360.} , {2,4,6,6,4,2} , false ),
    //      grid::StructuredGrid::YSpace( grid::LinearSpacing( {80.,-80.}, 6 ) ),
    //      Projection(),
    //      Domain() );

    std::string gridname = eckit::Resource<std::string>("--grid", "O8");

    StructuredGrid grid(gridname);

    int halo = eckit::Resource<int>("--halo", 2);
    util::Config config;
    config.set("halo", halo);
    config.set("periodic_points", true);
    functionspace::StructuredColumns fs(grid, grid::Partitioner("equal_regions"), config);

    Field field = fs.createField<long>(option::name("field"));

    auto value = array::make_view<long, 1>(field);
    auto xy    = array::make_view<double, 2>(fs.xy());
    auto g     = array::make_view<gidx_t, 1>(fs.global_index());
    auto r     = array::make_view<idx_t, 1>(fs.remote_index());
    auto p     = array::make_view<int, 1>(fs.partition());

    for (idx_t j = fs.j_begin(); j < fs.j_end(); ++j) {
        for (idx_t i = fs.i_begin(j); i < fs.i_end(j); ++i) {
            idx_t n  = fs.index(i, j);
            value(n) = util::microdeg(xy(n, XX));
        }
    }

    // EXPECT( fs.checksum(field) == "cef2694016492d408fa157b7c59ce741" );

    fs.haloExchange(field);

    // EXPECT( fs.checksum(field) == "cef2694016492d408fa157b7c59ce741" );

    ATLAS_TRACE_SCOPE("Output python") {
        eckit::PathName filepath("test_functionspace_StructuredColumns_halo_p" + std::to_string(mpi::comm().rank()) +
                                 ".py");

        std::ofstream f(filepath.asString().c_str(), std::ios::trunc);

        f << "\n"
             "import matplotlib.pyplot as plt"
             "\n"
             "from matplotlib.path import Path"
             "\n"
             "import matplotlib.patches as patches"
             "\n"
             ""
             "\n"
             "from itertools import cycle"
             "\n"
             "import matplotlib.cm as cm"
             "\n"
             "import numpy as np"
             "\n"
             ""
             "\n"
             "fig = plt.figure(figsize=(20,10))"
             "\n"
             "ax = fig.add_subplot(111,aspect='equal')"
             "\n"
             "";

        double xmin = std::numeric_limits<double>::max();
        double xmax = -std::numeric_limits<double>::max();
        double ymin = std::numeric_limits<double>::max();
        double ymax = -std::numeric_limits<double>::max();
        f << "\n"
             "x = [";
        for (idx_t j = fs.j_begin_halo(); j < fs.j_end_halo(); ++j) {
            for (idx_t i = fs.i_begin_halo(j); i < fs.i_end_halo(j); ++i) {
                idx_t n = fs.index(i, j);
                f << xy(n, XX) << ", ";
                xmin = std::min(xmin, xy(n, XX));
                xmax = std::max(xmax, xy(n, XX));
            }
        }
        f << "]";

        f << "\n"
             "y = [";
        for (idx_t j = fs.j_begin_halo(); j < fs.j_end_halo(); ++j) {
            for (idx_t i = fs.i_begin_halo(j); i < fs.i_end_halo(j); ++i) {
                idx_t n = fs.index(i, j);
                f << xy(n, YY) << ", ";
                ymin = std::min(ymin, xy(n, YY));
                ymax = std::max(ymax, xy(n, YY));
            }
        }
        f << "]";

        f << "\n"
             "g = [";
        for (idx_t j = fs.j_begin_halo(); j < fs.j_end_halo(); ++j) {
            for (idx_t i = fs.i_begin_halo(j); i < fs.i_end_halo(j); ++i) {
                idx_t n = fs.index(i, j);
                f << g(n) << ", ";
            }
        }
        f << "]";

        f << "\n"
             "p = [";
        for (idx_t j = fs.j_begin_halo(); j < fs.j_end_halo(); ++j) {
            for (idx_t i = fs.i_begin_halo(j); i < fs.i_end_halo(j); ++i) {
                idx_t n = fs.index(i, j);
                f << p(n) << ", ";
            }
        }
        f << "]";

        f << "\n"
             "r = [";
        for (idx_t j = fs.j_begin_halo(); j < fs.j_end_halo(); ++j) {
            for (idx_t i = fs.i_begin_halo(j); i < fs.i_end_halo(j); ++i) {
                idx_t n = fs.index(i, j);
                f << r(n) << ", ";
            }
        }
        f << "]";

        f << "\n"
             ""
             "\n"
             "c = [ cm.Paired( float(pp%13)/12. ) for pp in p ]"
             "\n"
             "ax.scatter(x, y, color=c, marker='o')"
             "\n"
             "for i in range("
          << fs.size()
          << "):"
             "\n"
             "  ax.annotate(g[i], (x[i],y[i]), fontsize=8)"
             "\n"
             "";
        f << "\n"
             "ax.set_xlim( "
          << std::min(0., xmin) << "-5, " << std::max(360., xmax)
          << "+5)"
             "\n"
             "ax.set_ylim( "
          << std::min(-90., ymin) << "-5, " << std::max(90., ymax)
          << "+5)"
             "\n"
             "ax.set_xticks([0,45,90,135,180,225,270,315,360])"
             "\n"
             "ax.set_yticks([-90,-45,0,45,90])"
             "\n"
             "plt.grid()"
             "\n"
             "plt.show()"
             "\n";
    }
}

CASE("StructuredColumns checksum matches BlockStructuredColumns for identical data") {
    using Value = double;
    mpi::Scope scope("self");

    auto grid = StructuredGrid("O8");
    auto structured_fs = [&] {
        util::Config config;
        config.set("halo", 0);
        return functionspace::StructuredColumns(grid, config);
    }();
    auto block_fs = [&] {
        util::Config config;
        config.set("halo", 0);
        config.set("nproma", 5);
        return functionspace::BlockStructuredColumns(grid, config);
    }();

    auto check_field_checksum_stable_with_blocked = [&](const util::Config& options, Value first_value) {
        Field structured = structured_fs.createField<Value>(option::name("structured") | options);
        Field blocked = block_fs.createField<Value>(option::name("blocked") | options);
        fill_structured_field(structured, first_value);
        block_fs.scatter(structured, blocked);
        EXPECT_EQ(structured_fs.checksum(structured), block_fs.checksum(blocked));
    };

    check_field_checksum_stable_with_blocked(option::variables(3) | option::levels(4), Value{1});
    check_field_checksum_stable_with_blocked(option::variables(3), Value{1000});
    check_field_checksum_stable_with_blocked(option::levels(4), Value{2000});
    check_field_checksum_stable_with_blocked(util::Config{}, Value{3000});

    auto check_fieldset_checksum_stable_with_blocked = [&](const util::Config& options_a, const util::Config& options_b) {
        FieldSet structured_fieldset;
        FieldSet blocked_fieldset;
        Field structured_a = structured_fs.createField<Value>(option::name("structured_a") | options_a);
        Field structured_b = structured_fs.createField<Value>(option::name("structured_b") | options_b);
        Field blocked_a = block_fs.createField<Value>(option::name("blocked_a") | options_a);
        Field blocked_b = block_fs.createField<Value>(option::name("blocked_b") | options_b);

        fill_structured_field(structured_a, Value{4000});
        fill_structured_field(structured_b, Value{5000});
        block_fs.scatter(structured_a, blocked_a);
        block_fs.scatter(structured_b, blocked_b);

        structured_fieldset.add(structured_a);
        structured_fieldset.add(structured_b);
        blocked_fieldset.add(blocked_a);
        blocked_fieldset.add(blocked_b);

        EXPECT_EQ(structured_fs.checksum(structured_fieldset), block_fs.checksum(blocked_fieldset));
    };
    check_fieldset_checksum_stable_with_blocked(/*a*/ option::variables(2) | option::levels(3), /*b*/ option::variables(4));
}

CASE("StructuredColumns checksum matches BlockStructuredColumns in MPI for identical scattered data") {
    using Value = double;

    constexpr int root = 0;
    auto grid = StructuredGrid("O8");
    auto structured_fs = [&] {
        util::Config config;
        config.set("halo", 0);
        return functionspace::StructuredColumns(grid, config);
    }();
    auto block_fs = [&] {
        util::Config config;
        config.set("halo", 0);
        config.set("nproma", 5);
        return functionspace::BlockStructuredColumns(grid, config);
    }();

    auto run_case = [&](const util::Config& options, Value first_value) {
        Field global = structured_fs.createField<Value>(option::name("global") | option::global(root) | options);
        Field structured = structured_fs.createField<Value>(option::name("structured") | options);
        Field blocked = block_fs.createField<Value>(option::name("blocked") | options);
        if (mpi::comm().rank() == root) {
            fill_structured_field(global, first_value);
        }
        structured_fs.scatter(global, structured);
        block_fs.scatter(global, blocked);
        EXPECT_EQ(structured_fs.checksum(structured), block_fs.checksum(blocked));
    };

    run_case(option::variables(3) | option::levels(4), Value{1});
    run_case(option::variables(3), Value{1000});
    run_case(option::levels(4), Value{2000});
    run_case(util::Config{}, Value{3000});

    Field global_a = structured_fs.createField<Value>(option::name("global_a") | option::global(root) |
                                                      option::variables(2) | option::levels(3));
    Field global_b = structured_fs.createField<Value>(option::name("global_b") | option::global(root) |
                                                      option::variables(4));
    Field structured_a = structured_fs.createField<Value>(option::name("structured_a") | option::variables(2) |
                                                          option::levels(3));
    Field structured_b = structured_fs.createField<Value>(option::name("structured_b") | option::variables(4));
    Field blocked_a = block_fs.createField<Value>(option::name("blocked_a") | option::variables(2) | option::levels(3));
    Field blocked_b = block_fs.createField<Value>(option::name("blocked_b") | option::variables(4));

    if (mpi::comm().rank() == root) {
        fill_structured_field(global_a, Value{4000});
        fill_structured_field(global_b, Value{5000});
    }

    structured_fs.scatter(global_a, structured_a);
    structured_fs.scatter(global_b, structured_b);
    block_fs.scatter(global_a, blocked_a);
    block_fs.scatter(global_b, blocked_b);

    FieldSet structured_fieldset;
    FieldSet blocked_fieldset;
    structured_fieldset.add(structured_a);
    structured_fieldset.add(structured_b);
    blocked_fieldset.add(blocked_a);
    blocked_fieldset.add(blocked_b);

    EXPECT_EQ(structured_fs.checksum(structured_fieldset), block_fs.checksum(blocked_fieldset));
}

CASE("StructuredColumns and BlockStructuredColumns checksums are stable with halo=2") {
    using Value = double;
    mpi::Scope scope("self");

    auto grid = StructuredGrid("O8");
    auto structured_fs = [&] {
        util::Config config;
        config.set("halo", 2);
        return functionspace::StructuredColumns(grid, config);
    }();
    auto block_fs = [&] {
        util::Config config;
        config.set("halo", 2);
        config.set("nproma", 5);
        return functionspace::BlockStructuredColumns(grid, config);
    }();
    util::Config options = option::variables(3) | option::levels(4);
    Field structured = structured_fs.createField<Value>(option::name("structured") | options);
    fill_structured_field(structured, Value{6000});
    FieldSet structured_fieldset;
    structured_fieldset.add(structured);
    EXPECT_EQ(structured_fs.checksum(structured), structured_fs.checksum(structured_fieldset));
    EXPECT_EQ(structured_fs.checksum(structured), structured_fs.checksum(structured));

    Field blocked = block_fs.createField<Value>(option::name("blocked") | options);
    block_fs.scatter(structured, blocked);
    FieldSet blocked_fieldset;
    blocked_fieldset.add(blocked);
    EXPECT_EQ(block_fs.checksum(blocked), block_fs.checksum(blocked_fieldset));
    EXPECT_EQ(block_fs.checksum(blocked), block_fs.checksum(blocked));
}

CASE("StructuredColumns checksum with halo=2 matches halo=0 exactly (mpi-serial)") {
    using Value = double;
    mpi::Scope scope("self");

    auto grid = StructuredGrid("O8");
    auto fs_h0 = [&] {
        util::Config config;
        config.set("halo", 0);
        return functionspace::StructuredColumns(grid, config);
    }();
    auto fs_h2 = [&] {
        util::Config config;
        config.set("halo", 2);
        return functionspace::StructuredColumns(grid, config);
    }();

    auto check_field_checksum_stable_with_halo = [&](const util::Config& options, Value first_value) {
        Field field_h0 = fs_h0.createField<Value>(option::name("field_h0") | options);
        Field field_h2 = fs_h2.createField<Value>(option::name("field_h2") | options);
        fill_structured_field(field_h0, first_value);
        fill_structured_field(field_h2, first_value);
        EXPECT_EQ(fs_h0.checksum(field_h0), fs_h2.checksum(field_h2));
    };

    check_field_checksum_stable_with_halo(option::variables(3) | option::levels(4), Value{12000});
    check_field_checksum_stable_with_halo(option::variables(3), Value{13000});
    check_field_checksum_stable_with_halo(option::levels(4), Value{14000});
    check_field_checksum_stable_with_halo(util::Config{}, Value{15000});

    auto check_fieldset_checksum_stable_with_halo = [&](const util::Config& options_a, const util::Config& options_b) {
        Field field_h0_a = fs_h0.createField<Value>(option::name("field_h0_a") | options_a);
        Field field_h0_b = fs_h0.createField<Value>(option::name("field_h0_b") | options_b);
        Field field_h2_a = fs_h2.createField<Value>(option::name("field_h2_a") | options_a);
        Field field_h2_b = fs_h2.createField<Value>(option::name("field_h2_b") | options_b);

        fill_structured_field(field_h0_a, Value{16000});
        fill_structured_field(field_h0_b, Value{17000});
        fill_structured_field(field_h2_a, Value{16000});
        fill_structured_field(field_h2_b, Value{17000});

        FieldSet fieldset_h0;
        FieldSet fieldset_h2;
        fieldset_h0.add(field_h0_a);
        fieldset_h0.add(field_h0_b);
        fieldset_h2.add(field_h2_a);
        fieldset_h2.add(field_h2_b);

        EXPECT_EQ(fs_h0.checksum(fieldset_h0), fs_h2.checksum(fieldset_h2));
    };

    check_fieldset_checksum_stable_with_halo(/*a*/ option::variables(2) | option::levels(3), /*b*/ option::variables(4));

}

CASE("StructuredColumns checksum with halo=2 matches halo=0 exactly (MPI)") {
    using Value = double;

    constexpr int root = 0;
    auto grid = StructuredGrid("O8");
    auto fs_h0 = [&] {
        return functionspace::StructuredColumns(grid, option::halo(0));
    }();
    auto fs_h2 = [&] {
        return functionspace::StructuredColumns(grid, option::halo(2));
    }();

    auto check_field_checksum_stable_with_halo = [&](const util::Config& options, Value first_value) {
        Field global = fs_h0.createField<Value>(option::name("global") | option::global(root) | options);
        Field local_h0 = fs_h0.createField<Value>(option::name("local_h0") | options);
        Field local_h2 = fs_h2.createField<Value>(option::name("local_h2") | options);
        if (mpi::comm().rank() == root) {
            fill_structured_field(global, first_value);
        }
        fs_h0.scatter(global, local_h0);
        fs_h2.scatter(global, local_h2);
        EXPECT_EQ(fs_h0.checksum(local_h0), fs_h2.checksum(local_h2));
    };

    check_field_checksum_stable_with_halo(option::variables(3) | option::levels(4), Value{18000});
    check_field_checksum_stable_with_halo(option::variables(3), Value{19000});
    check_field_checksum_stable_with_halo(option::levels(4), Value{20000});
    check_field_checksum_stable_with_halo(util::Config{}, Value{21000});

    auto check_fieldset_checksum_stable_with_halo = [&](const util::Config& options_a, const util::Config& options_b) {
        Field global_a = fs_h0.createField<Value>(option::name("global_a") | option::global(root) | options_a);
        Field global_b = fs_h0.createField<Value>(option::name("global_b") | option::global(root) | options_b);
        Field local_h0_a = fs_h0.createField<Value>(option::name("local_h0_a") | options_a);
        Field local_h0_b = fs_h0.createField<Value>(option::name("local_h0_b") | options_b);
        Field local_h2_a = fs_h2.createField<Value>(option::name("local_h2_a") | options_a);
        Field local_h2_b = fs_h2.createField<Value>(option::name("local_h2_b") | options_b);

        if (mpi::comm().rank() == root) {
            fill_structured_field(global_a, Value{22000});
            fill_structured_field(global_b, Value{23000});
        }

        fs_h0.scatter(global_a, local_h0_a);
        fs_h0.scatter(global_b, local_h0_b);
        fs_h2.scatter(global_a, local_h2_a);
        fs_h2.scatter(global_b, local_h2_b);

        FieldSet fieldset_h0;
        FieldSet fieldset_h2;
        fieldset_h0.add(local_h0_a);
        fieldset_h0.add(local_h0_b);
        fieldset_h2.add(local_h2_a);
        fieldset_h2.add(local_h2_b);

        EXPECT_EQ(fs_h0.checksum(fieldset_h0), fs_h2.checksum(fieldset_h2));
    };
    check_fieldset_checksum_stable_with_halo(/*a*/ option::variables(2) | option::levels(3), /*b*/ option::variables(4));
}

//-----------------------------------------------------------------------------

CASE("test_functionspace_StructuredColumns_halo checks without output") {
    std::string gridname = eckit::Resource<std::string>("--grid", "O8");

    StructuredGrid grid(gridname);

    int halo = eckit::Resource<int>("--halo", 2);
    util::Config config;
    config.set("halo", halo);
    config.set("levels", 10);
    config.set("periodic_points", true);
    functionspace::StructuredColumns fs(grid, grid::Partitioner("equal_regions"), config);

    Field field = fs.createField<long>(option::name("field"));

    auto value = array::make_view<long, 2>(field);
    auto xy    = array::make_view<double, 2>(fs.xy());

    for (idx_t j = fs.j_begin(); j < fs.j_end(); ++j) {
        for (idx_t i = fs.i_begin(j); i < fs.i_end(j); ++i) {
            idx_t n = fs.index(i, j);
            for (idx_t k = 0; k < fs.levels(); ++k) {
                value(n, k) = util::microdeg(xy(n, XX));
            }
        }
    }

    ATLAS_TRACE_SCOPE("control each value ")
    fs.parallel_for([&](idx_t n, idx_t k) { EXPECT(value(n, k) == util::microdeg(xy(n, XX))); });
    fs.parallel_for([&](idx_t n, idx_t i, idx_t j, idx_t k) { EXPECT(value(n, k) == util::microdeg(grid.x(i, j))); });

    Field fieldg = fs.createField(field, option::global());
    fs.gather(field, fieldg);

    ATLAS_TRACE_SCOPE("control_global") {
        auto valueg = array::make_view<long, 2>(fieldg);
        fs.parallel_for(option::global(), [&](idx_t n, idx_t i, idx_t j, idx_t k) {
            EXPECT(valueg(n, k) == util::microdeg(grid.x(i, j)));
        });
    }
}

//-----------------------------------------------------------------------------

CASE("test_functionspace_StructuredColumns halo exchange registration") {
    // Test by observing log with ATLAS_DEBUG=1
    // The HaloExchange Cache should be created twice, already found 4 times,
    // erased twice.

    std::string gridname = eckit::Resource<std::string>("--grid", "O8");

    StructuredGrid grid(gridname);

    util::Config config;
    config.set("levels", 10);
    config.set("periodic_points", true);
    for (idx_t i = 0; i < 3; ++i) {
        config.set("halo", 2);
        functionspace::StructuredColumns fs1(grid, grid::Partitioner("equal_regions"), config);
        config.set("halo", 4);
        functionspace::StructuredColumns fs2(grid, grid::Partitioner("equal_regions"), config);

        Field field1 = fs1.createField<long>(option::name("field"));
        Field field2 = fs2.createField<long>(option::name("field"));

        field1.haloExchange();
        field2.haloExchange();
    }
}


//-----------------------------------------------------------------------------

long innerproductwithhalo(const atlas::Field& f1, const atlas::Field& f2) {
    long sum(0);

    auto view1 = atlas::array::make_view<long, 2>(f1);
    auto view2 = atlas::array::make_view<long, 2>(f2);

    for (atlas::idx_t jn = 0; jn < f1.shape(0); ++jn) {
        for (atlas::idx_t jl = 0; jl < f1.levels(); ++jl) {
            sum += view1(jn, jl) * view2(jn, jl);
        }
    }

    atlas::mpi::comm().allReduceInPlace(sum, eckit::mpi::sum());
    return sum;
}


CASE("test_functionspace_StructuredColumns halo exchange adjoint test 1") {
    // Adjoint test for fields

    std::string gridname = eckit::Resource<std::string>("--grid", "S20x3");

    StructuredGrid grid(gridname);

    util::Config config;
    config.set("levels", 1);
    config.set("halo", 1);

    long sum1(0), sum2(0);
    functionspace::StructuredColumns fs(grid, grid::Partitioner("checkerboard"), config);
    Field fieldInit = fs.createField<long>(option::name("fieldInit"));

    // Field setup  values 1 in interior and zeros in halos.
    auto view1   = atlas::array::make_view<long, 2>(fieldInit);
    auto i_index = atlas::array::make_view<atlas::idx_t, 1>(fs.index_i());
    auto j_index = atlas::array::make_view<atlas::idx_t, 1>(fs.index_j());
    for (atlas::idx_t jn = 0; jn < fs.sizeOwned(); ++jn) {
        for (atlas::idx_t jl = 0; jl < fieldInit.levels(); ++jl) {
            view1(jn, jl) = 1;
            std::cout << " initial  interior regions:: "
                      << " size = " << atlas::mpi::comm().size() << " "
                      << " jn = " << jn << " "
                      << " rank = " << atlas::mpi::comm().rank() << " "
                      << " jn = " << jn << " "
                      << " view1 = " << view1(jn, 0) << " "
                      << " i_index = " << i_index(jn) << " "
                      << " j_index = " << j_index(jn) << std::endl;
        }
    }

    Field fieldTmp = fs.createField<long>(option::name("fieldTmp"));
    auto viewTmp   = atlas::array::make_view<long, 2>(fieldTmp);
    for (atlas::idx_t jn = 0; jn < fieldTmp.shape(0); ++jn) {
        for (atlas::idx_t jl = 0; jl < fieldTmp.levels(); ++jl) {
            viewTmp(jn, jl) = view1(jn, jl);
        }
    }
    fieldTmp.haloExchange();

    sum1 = innerproductwithhalo(fieldTmp, fieldTmp);

    fieldTmp.adjointHaloExchange();

    sum2 = innerproductwithhalo(fieldInit, fieldTmp);

    atlas::Log::info() << "sum1 " << sum1 << "sum2 " << sum2 << std::endl;
    EXPECT(sum1 == sum2);

    // check that it is zero is in the halo
    sum1 = 0;
    for (atlas::idx_t jn = fs.sizeOwned(); jn < fieldTmp.shape(0); ++jn) {
        for (atlas::idx_t jl = 0; jl < fieldTmp.levels(); ++jl) {
            sum1 += viewTmp(jn, jl);
        }
    }
    atlas::mpi::comm().allReduceInPlace(sum1, eckit::mpi::sum());
    EXPECT(sum1 == 0);
}


CASE("create_aligned_field") {
    std::string gridname = eckit::Resource<std::string>("--grid", "S20x3");
    Grid grid(gridname);
    functionspace::StructuredColumns fs(grid, option::levels(5));
    Field field      = fs.createField<double>(option::variables(3) | option::alignment(4));
    auto check_field = [&](const Field& field) {
        EXPECT_EQ(field.shape()[0], fs.size());
        EXPECT_EQ(field.shape()[1], 5);
        EXPECT_EQ(field.shape()[2], 3);
        EXPECT_EQ(field.size(), fs.size() * 5 * 3);
        EXPECT_EQ(field.contiguous(), false);
        EXPECT_EQ(field.strides()[0], 5 * 4);
        EXPECT_EQ(field.strides()[1], 4);
        EXPECT_EQ(field.strides()[2], 1);
    };
    check_field(field);
}

//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
