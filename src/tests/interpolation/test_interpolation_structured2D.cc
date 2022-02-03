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

using atlas::functionspace::NodeColumns;
using atlas::functionspace::StructuredColumns;
using atlas::util::Config;

namespace atlas {
namespace test {

//-----------------------------------------------------------------------------

static Config scheme() {
    Config scheme;
    std::string scheme_str = eckit::Resource<std::string>("--scheme", "linear");
    if (scheme_str == "linear") {
        scheme.set("type", "structured-linear2D");
        scheme.set("halo", 1);
        // The stencil does not require any halo, but we set it to 1 for pole treatment!
    }
    if (scheme_str == "cubic") {
        scheme.set("type", "structured-cubic2D");
        scheme.set("halo", 2);
    }
    if (scheme_str == "quasicubic") {
        scheme.set("type", "structured-quasicubic2D");
        scheme.set("halo", 2);
    }
    scheme.set("name", scheme_str);
    return scheme;
}

std::string input_gridname(const std::string& default_grid) {
    return eckit::Resource<std::string>("--input-grid", default_grid);
}

std::string output_gridname(const std::string& default_grid) {
    return eckit::Resource<std::string>("--output-grid", default_grid);
}

Grid rotated_mercator() {
    Config gridspec;
    gridspec.set("type", "regional");
    gridspec.set("nx", 50);
    gridspec.set("ny", 40);
    gridspec.set("dx", 50000);
    gridspec.set("dy", 50000);
    gridspec.set("y_numbering", -1);
    gridspec.set("lonlat(centre)", std::vector<double>{4., 50});
    gridspec.set("projection", [] {
        Config projection;
        projection.set("type", "rotated_mercator");
        projection.set("north_pole", std::vector<double>{-176., 40.});
        return projection;
    }());
    return Grid{gridspec};
}

Grid lambert() {
    Config gridspec;
    gridspec.set("type", "regional");
    gridspec.set("nx", 50);
    gridspec.set("ny", 40);
    gridspec.set("dx", 50000);
    gridspec.set("dy", 50000);
    gridspec.set("y_numbering", -1);
    gridspec.set("lonlat(centre)", std::vector<double>{4., 50});
    gridspec.set("projection", [] {
        Config projection;
        projection.set("type", "lambert_conformal_conic");
        projection.set("longitude0", 4.);
        projection.set("latitude0", 50.);
        return projection;
    }());
    return Grid{gridspec};
}

Grid rotated(const std::string& name) {
    Config gridspec;
    gridspec.set("name", name);
    gridspec.set("projection", [] {
        Config projection;
        projection.set("type", "rotated_lonlat");
        projection.set("north_pole", std::vector<double>{-176., 40.});
        return projection;
    }());
    return Grid{gridspec};
}

FunctionSpace output_functionspace(const Grid& grid) {
    MeshGenerator meshgen("structured");
    Mesh output_mesh = meshgen.generate(grid);
    return NodeColumns{output_mesh};
}

CASE("which scheme?") {
    Log::info() << scheme().getString("type") << std::endl;
}

CASE("test_interpolation_structured using functionspace API") {
    Grid grid(input_gridname("O32"));

    // Cubic interpolation requires a StructuredColumns functionspace with 2 halos
    StructuredColumns input_fs(grid, scheme());

    auto test = [&](const FunctionSpace& output_fs) {
        // The output functionspace can currently be either NodeColumns or PointCloud

        Interpolation interpolation(scheme(), input_fs, output_fs);

        Field field_source = input_fs.createField<double>(option::name("source"));
        Field field_target = output_fs.createField<double>(option::name("target"));

        auto lonlat = array::make_view<double, 2>(input_fs.xy());
        auto source = array::make_view<double, 1>(field_source);
        for (idx_t n = 0; n < input_fs.size(); ++n) {
            source(n) = util::function::vortex_rollup(lonlat(n, LON), lonlat(n, LAT), 1.);
        }

        EXPECT(field_source.dirty());

        interpolation.execute(field_source, field_target);
    };

    SECTION("Interpolate from " + grid.name() + " to " + output_gridname("O64")) {
        EXPECT_NO_THROW(test(output_functionspace(Grid{output_gridname("O64")})));
    }
    SECTION("Interpolate from " + grid.name() + " to rotated " + output_gridname("O64")) {
        EXPECT_NO_THROW(test(output_functionspace(rotated(output_gridname("O64")))));
    }
    SECTION("Interpolate from " + grid.name() + " to lambert") {
        EXPECT_NO_THROW(test(output_functionspace(lambert())));
    }
    SECTION("Interpolate from " + grid.name() + " to rotaded_mercator") {
        EXPECT_NO_THROW(test(output_functionspace(rotated_mercator())));
    }
}

CASE("test_interpolation_structured using grid API") {
    // Using the grid API we can hide interpolation method specific requirements
    // such as which functionspace needs to be set-up.
    // Currently the assumption is that grids are serial

    Grid input_grid(input_gridname("O32"));

    auto test = [&](const Grid& output_grid) {
        Interpolation interpolation(scheme(), input_grid, output_grid);

        // Allocate and initialise own memory here to show possibilities
        // Note that allocated size must be possibly enlarged depending on interpolation method
        // allocates stencil halo
        std::vector<double> src_data(interpolation.source().size());
        std::vector<double> tgt_data(interpolation.target().size());

        idx_t n{0};
        for (auto p : input_grid.lonlat()) {
            src_data[n++] = util::function::vortex_rollup(p.lon(), p.lat(), 1.);
        }

        // Wrap memory in atlas Fields and interpolate
        Field field_source{"source", src_data.data(), array::make_shape(src_data.size())};
        Field field_target{"target", tgt_data.data(), array::make_shape(tgt_data.size())};

        // Wrapping data does not set the field to have dirty halo's
        {
            EXPECT(not field_source.dirty());
            field_source.set_dirty();
        }

        EXPECT(field_source.dirty());
        interpolation.execute(field_source, field_target);

        ATLAS_TRACE_SCOPE("output") {
            output::Gmsh gmsh(scheme().getString("name") + "-output-section" + std::to_string(_subsection) + ".msh",
                              Config("coordinates", "xy"));
            gmsh.write(MeshGenerator("structured").generate(output_grid));
            gmsh.write(field_target, StructuredColumns(output_grid));
        }
    };

    SECTION("Interpolate from " + input_grid.name() + " to " + output_gridname("O64")) {
        EXPECT_NO_THROW(test(Grid{output_gridname("O64")}));
    }
    SECTION("Interpolate from " + input_grid.name() + " to rotated " + output_gridname("O64")) {
        EXPECT_NO_THROW(test(rotated(output_gridname("O64"))));
    }
    SECTION("Interpolate from " + input_grid.name() + " to lambert") { EXPECT_NO_THROW(test(lambert())); }
    SECTION("Interpolate from " + input_grid.name() + " to rotaded_mercator") {
        EXPECT_NO_THROW(test(rotated_mercator()));
    }
}


CASE("test_interpolation_structured using fs API multiple levels") {
    Grid input_grid(input_gridname("O32"));
    Grid output_grid(output_gridname("O64"));

    // Cubic interpolation requires a StructuredColumns functionspace with 2 halos
    StructuredColumns input_fs(input_grid, scheme() | option::levels(3));

    MeshGenerator meshgen("structured");
    Mesh output_mesh        = meshgen.generate(output_grid);
    FunctionSpace output_fs = NodeColumns{output_mesh, option::levels(3)};

    Interpolation interpolation(scheme(), input_fs, output_fs);

    Field field_source = input_fs.createField<double>(option::name("source"));
    Field field_target = output_fs.createField<double>(option::name("target"));

    auto lonlat = array::make_view<double, 2>(input_fs.xy());
    auto source = array::make_view<double, 2>(field_source);
    for (idx_t n = 0; n < input_fs.size(); ++n) {
        for (idx_t k = 0; k < 3; ++k) {
            source(n, k) = util::function::vortex_rollup(lonlat(n, LON), lonlat(n, LAT), 0.5 + double(k) / 2);
        }
    }
    interpolation.execute(field_source, field_target);

    ATLAS_TRACE_SCOPE("output") {
        output::Gmsh gmsh(
            scheme().getString("name") + "-multilevel-output-section" + std::to_string(_subsection) + ".msh",
            Config("coordinates", "xy"));
        gmsh.write(output_mesh);
        output_fs.haloExchange(field_target);
        gmsh.write(field_target);
    }
}

template <typename Value>
struct AdjointTolerance {
    static const Value value;
};
template <>
const double AdjointTolerance<double>::value = 2.e-14;
template <>
const float AdjointTolerance<float>::value = 2.e-5;


template <typename Value>
void test_interpolation_structured_using_fs_API_for_fieldset() {
    Grid input_grid(input_gridname("O32"));
    Grid output_grid(output_gridname("O64"));

    // Cubic interpolation requires a StructuredColumns functionspace with 2 halos
    StructuredColumns input_fs(input_grid, scheme() | option::levels(3));

    MeshGenerator meshgen("structured");
    Mesh output_mesh        = meshgen.generate(output_grid);
    FunctionSpace output_fs = NodeColumns{output_mesh, option::levels(3)};

    auto lonlat = array::make_view<double, 2>(input_fs.xy());

    FieldSet fields_source;
    FieldSet fields_target;
    for (idx_t f = 0; f < 3; ++f) {
        auto field_source = fields_source.add(input_fs.createField<Value>(option::name("field " + std::to_string(f))));
        fields_target.add(output_fs.createField<Value>(option::name("field " + std::to_string(f))));

        auto source = array::make_view<Value, 2>(field_source);
        for (idx_t n = 0; n < input_fs.size(); ++n) {
            for (idx_t k = 0; k < 3; ++k) {
                source(n, k) = util::function::vortex_rollup(lonlat(n, LON), lonlat(n, LAT), 0.5 + double(k) / 2);
            }
        }
    }

    ATLAS_TRACE_SCOPE("with matrix") {
        Interpolation interpolation(scheme(), input_fs, output_fs);
        interpolation.execute(fields_source, fields_target);

        ATLAS_TRACE_SCOPE("output") {
            output::Gmsh gmsh(scheme().getString("name") + "-multilevel-fieldset-output-with-matrix-" +
                                  array::make_datatype<Value>().str() + ".msh",
                              Config("coordinates", "xy"));
            gmsh.write(output_mesh);
            output_fs.haloExchange(fields_target);
            gmsh.write(fields_target);
        }
    }

    ATLAS_TRACE_SCOPE("with matrix adjoint") {
        Interpolation interpolation(scheme() | Config("adjoint", true), input_fs, output_fs);

        std::vector<Value> AxAx(fields_source.field_names().size(), 0.);
        std::vector<Value> xAtAx(fields_source.field_names().size(), 0.);

        FieldSet fields_source_reference;
        for (atlas::Field& field : fields_source) {
            Field temp_field(field.name(), field.datatype().kind(), field.shape());
            temp_field.set_levels(field.levels());

            auto fieldInView  = array::make_view<Value, 2>(field);
            auto fieldOutView = array::make_view<Value, 2>(temp_field);

            for (atlas::idx_t jn = 0; jn < temp_field.shape(0); ++jn) {
                for (atlas::idx_t jl = 0; jl < temp_field.levels(); ++jl) {
                    fieldOutView(jn, jl) = fieldInView(jn, jl);
                }
            }
            fields_source_reference.add(temp_field);
        }

        interpolation.execute(fields_source, fields_target);

        std::size_t fIndx(0);
        auto source_names = fields_source.field_names();
        for (const std::string& s : fields_target.field_names()) {
            auto target = array::make_view<Value, 2>(fields_target[s]);
            auto source = array::make_view<Value, 2>(fields_source[source_names[fIndx]]);

            for (idx_t n = 0; n < input_fs.size(); ++n) {
                for (idx_t k = 0; k < 3; ++k) {
                    AxAx[fIndx] += source(n, k) * source(n, k);
                }
            }

            for (idx_t n = 0; n < output_fs.size(); ++n) {
                for (idx_t k = 0; k < 3; ++k) {
                    AxAx[fIndx] += target(n, k) * target(n, k);
                }
            }
            fIndx += 1;
        }

        interpolation.execute_adjoint(fields_source, fields_target);

        fIndx = 0;
        for (const std::string& s : fields_source.field_names()) {
            auto source_reference = array::make_view<Value, 2>(fields_source_reference[s]);
            auto source           = array::make_view<Value, 2>(fields_source[s]);

            for (idx_t n = 0; n < input_fs.size(); ++n) {
                for (idx_t k = 0; k < 3; ++k) {
                    xAtAx[fIndx] += source(n, k) * source_reference(n, k);
                }
            }
            fIndx += 1;
        }

        for (std::size_t t = 0; t < AxAx.size(); ++t) {
            Log::debug() << " Adjoint test t  = " << t << " (Ax).(Ax) = " << AxAx[t] << " x.(AtAx) = " << xAtAx[t]
                         << " std::abs( 1.0 - xAtAx[t]/AxAx[t] ) " << std::abs(1.0 - xAtAx[t] / AxAx[t])
                         << " AdjointTolerance<Value>::value " << AdjointTolerance<Value>::value << std::endl;

            EXPECT(std::abs(1.0 - xAtAx[t] / AxAx[t]) < AdjointTolerance<Value>::value);
        }
    }

    ATLAS_TRACE_SCOPE("matrix free") {
        Interpolation interpolation(scheme() | Config("matrix_free", true), input_fs, output_fs);
        interpolation.execute(fields_source, fields_target);
        ATLAS_TRACE_SCOPE("output") {
            output::Gmsh gmsh(scheme().getString("name") + "-multilevel-fieldset-output-section-matrix-free-" +
                                  array::make_datatype<Value>().str() + ".msh",
                              Config("coordinates", "xy"));
            gmsh.write(output_mesh);
            output_fs.haloExchange(fields_target);
            gmsh.write(fields_target);
        }
    }
}

CASE("test_interpolation_structured using fs API for fieldset (Value=double)") {
    test_interpolation_structured_using_fs_API_for_fieldset<double>();
}

CASE("test_interpolation_structured using fs API for fieldset (Value=float)") {
    test_interpolation_structured_using_fs_API_for_fieldset<float>();
}


/// @brief Compute magnitude of flow with rotation-angle beta
/// (beta=0 --> zonal, beta=pi/2 --> meridional)
Field rotated_flow(const StructuredColumns& fs, const double& beta) {
    const double radius  = util::Earth::radius();
    const double USCAL   = 20.;
    const double pvel    = USCAL / radius;
    const double deg2rad = M_PI / 180.;

    array::ArrayView<double, 2> lonlat_deg = array::make_view<double, 2>(fs.xy());

    Field field                     = fs.createField<double>(option::vector());
    array::ArrayView<double, 2> var = array::make_view<double, 2>(field);

    idx_t nnodes = var.shape(0);
    for (idx_t jnode = 0; jnode < nnodes; ++jnode) {
        double x        = lonlat_deg(jnode, LON) * deg2rad;
        double y        = lonlat_deg(jnode, LAT) * deg2rad;
        double Ux       = pvel * (std::cos(beta) + std::tan(y) * std::cos(x) * std::sin(beta)) * radius * std::cos(y);
        double Uy       = -pvel * std::sin(x) * std::sin(beta) * radius;
        var(jnode, LON) = Ux;
        var(jnode, LAT) = Uy;
    }
    return field;
}


CASE("test_interpolation_structured for vectors") {
    Grid input_grid(input_gridname("O32"));
    Grid output_grid(output_gridname("O64"));

    // Cubic interpolation requires a StructuredColumns functionspace with 2 halos
    StructuredColumns input_fs(input_grid, option::halo(2));

    MeshGenerator meshgen("structured");
    Mesh output_mesh        = meshgen.generate(output_grid);
    FunctionSpace output_fs = NodeColumns{output_mesh};

    Interpolation interpolation(scheme(), input_fs, output_fs);

    Field field_source = rotated_flow(input_fs, M_PI_4);
    Field field_target = output_fs.createField<double>(option::name("target") | option::vector());

    interpolation.execute(field_source, field_target);

    ATLAS_TRACE_SCOPE("output") {
        output::Gmsh gmsh(scheme().getString("name") + "-vector-output-section" + std::to_string(_subsection) + ".msh",
                          Config("coordinates", "xy"));
        gmsh.write(output_mesh);
        output_fs.haloExchange(field_target);
        gmsh.write(field_target);
    }
}

CASE("ATLAS-315: Target grid with domain West of 0 degrees Lon") {
    Grid input_grid(input_gridname("O256"));
    Grid output_grid(output_gridname("L45"), RectangularDomain({-10, 10}, {20, 60}));

    Interpolation interpolation(scheme() | Config("matrix_free", true), input_grid, output_grid);

    auto field_src = interpolation.source().createField<double>(option::name("src"));
    auto field_tgt = interpolation.target().createField<double>(option::name("tgt"));

    auto lonlat        = array::make_view<double, 2>(StructuredColumns(interpolation.source()).xy());
    auto source        = array::make_view<double, 1>(field_src);
    constexpr double k = 1;
    for (idx_t n = 0; n < source.size(); ++n) {
        source(n) = util::function::vortex_rollup(lonlat(n, LON), lonlat(n, LAT), 0.5 + double(k) / 2);
    }

    interpolation.execute(field_src, field_tgt);
}

}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
