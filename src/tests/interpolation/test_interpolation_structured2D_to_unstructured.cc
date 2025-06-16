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
#include "atlas/functionspace/NodeColumns.h"
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

FunctionSpace output_functionspace_match() {
    std::vector<PointXY> points;
    if (mpi::size() == 2) {
        if (mpi::rank() == 0) {
            points = std::vector<PointXY>{
                {45., 45.}, {90., 45.}, {135., 45.}, {180., 45.}, {225., 45.}, {270., 45.}, {315., 45.},
            };
        }
        if (mpi::rank() == 1) {
            points = std::vector<PointXY>{
                {45., -45.}, {90., -45.}, {135., -45.}, {180., -45.}, {225., -45.}, {270., -45.}, {315., -45.},
            };
        }
    }
    else if (mpi::size() == 1) {
        points = std::vector<PointXY>{
            {45., 45.},  {90., 45.},  {135., 45.},  {180., 45.},  {225., 45.},  {270., 45.},  {315., 45.},
            {45., -45.}, {90., -45.}, {135., -45.}, {180., -45.}, {225., -45.}, {270., -45.}, {315., -45.},
        };
    }
    else {
        ATLAS_NOTIMPLEMENTED;
    }
    return PointCloud(points);
}
FunctionSpace output_functionspace_nomatch() {
    std::vector<PointXY> points;
    if (mpi::size() == 2) {
        if (mpi::rank() == 0) {
            points = std::vector<PointXY>{
                {45., 45.}, {90., 45.}, {135., 45.}, {45., -45.}, {90., -45.}, {135., -45.},
            };
        }
        if (mpi::rank() == 1) {
            points = std::vector<PointXY>{
                {180., 45.},  {225., 45.},  {270., 45.},  {315., 45.},
                {180., -45.}, {225., -45.}, {270., -45.}, {315., -45.},
            };
        }
    }
    else if (mpi::size() == 1) {
        points = std::vector<PointXY>{
            {45., 45.},  {90., 45.},  {135., 45.},  {180., 45.},  {225., 45.},  {270., 45.},  {315., 45.},
            {45., -45.}, {90., -45.}, {135., -45.}, {180., -45.}, {225., -45.}, {270., -45.}, {315., -45.},
        };
    }
    else {
        ATLAS_NOTIMPLEMENTED;
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


CASE("test_match") {
    if(mpi::size() > 2) {
        return;
    }
    idx_t nb_fields = 2;
    idx_t nb_levels = 3;

    Grid input_grid(input_gridname("O32"));
    StructuredColumns input_fs(input_grid, option::halo(1) | option::levels(nb_levels));

    FunctionSpace output_fs = output_functionspace_match();

    Interpolation interpolation(option::type("structured-linear"), input_fs, output_fs);

    FieldSet fields_source = create_source_fields(input_fs, nb_fields, nb_levels);
    FieldSet fields_target = create_target_fields(output_fs, nb_fields, nb_levels);

    interpolation.execute(fields_source, fields_target);
}
CASE("test_nomatch") {
    if(mpi::size() > 2) {
        return;
    }

    idx_t nb_fields = 2;
    idx_t nb_levels = 3;

    Grid input_grid(input_gridname("O32"));
    StructuredColumns input_fs(input_grid, option::halo(1) | option::levels(nb_levels));

    FunctionSpace output_fs = output_functionspace_nomatch();

    if (false)  // expected to throw
    {
        Interpolation interpolation(option::type("structured-linear"), input_fs, output_fs);

        FieldSet fields_source = create_source_fields(input_fs, nb_fields, nb_levels);
        FieldSet fields_target = create_target_fields(output_fs, nb_fields, nb_levels);

        interpolation.execute(fields_source, fields_target);
    }
}

CASE("test_nomatch 2") {
    if(mpi::size() > 2) {
        return;
    }

    idx_t nb_fields = 2;
    idx_t nb_levels = 3;

    Grid input_grid(input_gridname("O32"));
    StructuredColumns input_fs(input_grid, option::halo(1) | option::levels(nb_levels));

    FunctionSpace output_fs = output_functionspace_nomatch();

    if (false)  // expected to throw
    {
        Interpolation interpolation(option::type("structured-linear"), input_fs, output_fs);

        FieldSet fields_source = create_source_fields(input_fs, nb_fields, nb_levels);
        FieldSet fields_target = create_target_fields(output_fs, nb_fields, nb_levels);

        interpolation.execute(fields_source, fields_target);
    }
}

}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
