/*
 * (C) Copyright 2026 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/array.h"
#include "atlas/field/FieldBuilder.h"
#include "atlas/field/Field.h"
#include "atlas/functionspace/StructuredColumns.h"
#include "atlas/grid/Grid.h"
#include "atlas/grid/StructuredGrid.h"
#include "atlas/util/Config.h"
#include "atlas/util/Geometry.h"

#include "tests/AtlasTestEnvironment.h"

namespace atlas {
namespace test {

namespace {

void validate_build_area_field(const Grid& grid) {
    const double one_square_metre = 1.e6;
    functionspace::StructuredColumns fs(grid, option::halo(1));

    field::FieldBuilder build_area_field{"grid-box-area", fs};
    Field area_field = build_area_field();
    auto area         = array::make_view<double, 1>(area_field);

    EXPECT_EQ(area_field.name(), "area");
    EXPECT_EQ(area_field.shape(0), fs.size());

    idx_t checked_rows = 0;
    for (idx_t j = fs.j_begin_halo(); j < fs.j_end_halo(); ++j) {
        const idx_t i_begin = fs.i_begin_halo(j);
        const idx_t i_end   = fs.i_end_halo(j);
        if (i_begin == i_end) {
            continue;
        }

        const double row_area = area(fs.index(i_begin, j));

        Log::info() << "Row " << j << " representative area = " << row_area / 1.e6 << " km^2" << std::endl;
        EXPECT(std::isfinite(row_area));
        EXPECT(row_area >= one_square_metre);
        ++checked_rows;

        for (idx_t i = i_begin; i < i_end; ++i) {
            const idx_t point_index = fs.index(i, j);
            EXPECT_EQ(area(point_index), row_area);
        }
    }
    EXPECT(checked_rows > 0);
}

Grid rotated_o32() {
    return Grid("O32", Projection(option::type("rotated_lonlat") |
                                  util::Config("north_pole", {-176., 40.})));
}

Grid regional_lonlat() {
    const int Nx = 8, Ny = 8;
    const double xmin = +20, xmax = +60, ymin = +20, ymax = +60;

    StructuredGrid::XSpace xspace(util::Config("type", "linear")("N", Nx)("start", xmin)("end", xmax));
    StructuredGrid::YSpace yspace(util::Config("type", "linear")("N", Ny)("start", ymin)("end", ymax));

    return StructuredGrid(xspace, yspace);
}

Grid regional_lambert() {
    util::Config gridspec;
    gridspec.set("type", "regional");
    gridspec.set("nx", 21);
    gridspec.set("ny", 31);
    gridspec.set("dx", 8000.);
    gridspec.set("dy", 9000.);
    gridspec.set("lonlat(centre)", {9.9, 56.3});
    gridspec.set("projection", [] {
        util::Config projection;
        projection.set("type", "lambert_conformal_conic");
        projection.set("longitude0", 0.);
        projection.set("latitude0", 56.3);
        return projection;
    }());
    return Grid(gridspec);
}


}  // namespace

CASE("test Grid(O32) area field") {
    validate_build_area_field(Grid("O32"));
}

CASE("test global regular lonlat 5 degree area field") {
    validate_build_area_field(Grid("L72x37"));
}

CASE("test rotated Grid(O32) area field") {
    validate_build_area_field(rotated_o32());
}

CASE("test regional lonlat area field") {
    validate_build_area_field(regional_lonlat());
}

CASE("test regional lambert area field") {
    validate_build_area_field(regional_lambert());
}

CASE("test grid-box-area configuration") {
    functionspace::StructuredColumns fs(Grid("O32"), option::halo(1));
    const idx_t point_index = fs.index(fs.i_begin(fs.j_begin()), fs.j_begin());

    Field default_area = field::FieldBuilder("grid-box-area", fs)();
    Field named_area = field::FieldBuilder("grid-box-area", fs, util::Config("field_name", "cell_area"))();
    Field radius_area = field::FieldBuilder("grid-box-area", fs, util::Config("radius", 2.))();
    Field unit_sphere_area = field::FieldBuilder("grid-box-area", fs, util::Config("geometry", "UnitSphere"))();
    Field geometry_overrides_radius = field::FieldBuilder("grid-box-area", fs, util::Config("radius", 2.)("geometry", "UnitSphere"))();

    const auto default_view = array::make_view<double, 1>(default_area);
    const auto radius_view = array::make_view<double, 1>(radius_area);
    const auto unit_sphere_view = array::make_view<double, 1>(unit_sphere_area);
    const auto override_view = array::make_view<double, 1>(geometry_overrides_radius);

    EXPECT_EQ(default_area.name(), "area");
    EXPECT_EQ(named_area.name(), "cell_area");

    const double earth_radius = geometry::Earth().radius();
    EXPECT_APPROX_EQ(radius_view(point_index) / default_view(point_index), 4. / (earth_radius * earth_radius), 1.e-14);
    EXPECT_APPROX_EQ(unit_sphere_view(point_index) / default_view(point_index),
                     1. / (earth_radius * earth_radius), 1.e-14);
    EXPECT_EQ(override_view(point_index), unit_sphere_view(point_index));
}

}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}