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
#include "atlas/field/Field.h"
#include "atlas/field/FieldSet.h"
#include "atlas/functionspace/PointCloud.h"
#include "atlas/functionspace/StructuredColumns.h"
#include "atlas/grid/Grid.h"
#include "atlas/grid/Iterator.h"
#include "atlas/interpolation.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/meshgenerator.h"
#include "atlas/output/Gmsh.h"
#include "atlas/util/CoordinateEnums.h"
#include "atlas/util/function/VortexRollup.h"
#include "atlas/util/PolygonXY.h"

#include "tests/AtlasTestEnvironment.h"

using atlas::functionspace::PointCloud;
using atlas::functionspace::StructuredColumns;
using atlas::util::Config;

namespace atlas {
namespace test {

//-----------------------------------------------------------------------------

std::string input_gridname(const std::string& default_grid) {
    return eckit::Resource<std::string>("--input-grid", default_grid);
}

FunctionSpace output_functionspace( const FunctionSpace& input_functionspace, bool expect_fail ) {

    std::vector<PointXY> all_points {
        {360., 90.},
        {360., 0.},
        {360., -90.},
        {0.,0.},
    };

    // Only keep points that match the input partitioning.
    // Note that with following implementation it could be that some points
    // are present in two partitions, but it is not a problem for this test purpose.
    std::vector<PointXY> points;
    auto polygon = util::PolygonXY{input_functionspace.polygon()};
    for (const auto& p : all_points ) {
        if (polygon.contains(p)) {
            points.emplace_back(p);
        }
    }
    if( expect_fail && mpi::rank() == mpi::size() - 1 ) {
        points.emplace_back(720,0.);
    }
    return PointCloud(points);
}


FieldSet create_source_fields(StructuredColumns& fs, idx_t nb_fields, idx_t nb_levels) {
    using Value = double;
    FieldSet fields_source;
    auto lonlat = array::make_view<double, 2>(fs.xy());
    for (idx_t f = 0; f < nb_fields; ++f) {
        auto field_source = fields_source.add(fs.createField<Value>());
        auto source       = array::make_view<Value, 2>(field_source);
        for (idx_t n = 0; n < fs.size(); ++n) {
            for (idx_t k = 0; k < nb_levels; ++k) {
                source(n, k) = util::function::vortex_rollup(lonlat(n, LON), lonlat(n, LAT), 0.5 + double(k) / 2);
            }
        };
    }
    return fields_source;
}
FieldSet create_target_fields(FunctionSpace& fs, idx_t nb_fields, idx_t nb_levels) {
    using Value = double;
    FieldSet fields_target;
    for (idx_t f = 0; f < nb_fields; ++f) {
        fields_target.add(fs.createField<Value>(option::levels(nb_levels)));
    }
    return fields_target;
}

void do_test( std::string type, int input_halo, bool matrix_free, bool expect_failure ) {
    idx_t nb_fields = 2;
    idx_t nb_levels = 3;

    Grid input_grid(input_gridname("O32"));
    StructuredColumns input_fs(input_grid, option::levels(nb_levels) |
        option::halo(input_halo));

    FunctionSpace output_fs = output_functionspace(input_fs, expect_failure);

    Interpolation interpolation(option::type(type) |
        util::Config("matrix_free",matrix_free) |
        util::Config("verbose",eckit::Resource<bool>("--verbose",false)),
        input_fs, output_fs);

    FieldSet fields_source = create_source_fields(input_fs, nb_fields, nb_levels);
    FieldSet fields_target = create_target_fields(output_fs, nb_fields, nb_levels);

    interpolation.execute(fields_source, fields_target);
}

CASE("test structured-bilinear, halo 2, with matrix") {
    EXPECT_NO_THROW( do_test("structured-linear2D",2,false,false) );
}

CASE("test structured-bilinear, halo 2, without matrix") {
    EXPECT_NO_THROW( do_test("structured-linear2D",2,true,false) );
}

CASE("test structured-bilinear, halo 2, with matrix, expected failure") {
    EXPECT_THROWS_AS( do_test("structured-linear2D",2,false,true), eckit::Exception );
}

CASE("test structured-bilinear, halo 2, without matrix, expected failure") {
    EXPECT_THROWS_AS( do_test("structured-linear2D",2,false,true), eckit::Exception );
}

CASE("test structured-bilinear, halo 1, with matrix, expected failure") {
    EXPECT_THROWS_AS( do_test("structured-linear2D",1,false,false), eckit::Exception );
}

CASE("test structured-bilinear, halo 1, without matrix, expected failure") {
    EXPECT_THROWS_AS( do_test("structured-linear2D",1,true,false), eckit::Exception );
}

CASE("test structured-bicubic, halo 3, with matrix") {
    EXPECT_NO_THROW( do_test("structured-cubic2D",3,false,false) );
}

CASE("test structured-bicubic, halo 2, with matrix") {
    EXPECT_THROWS_AS( do_test("structured-cubic2D",2,false,false), eckit::Exception );
}


}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
