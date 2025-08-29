/*
 * (C) Copyright 2024 ECMWF.
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
#include "atlas/functionspace/StructuredColumns.h"
#include "atlas/grid/Grid.h"
#include "atlas/grid/Iterator.h"
#include "atlas/interpolation.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/meshgenerator.h"
#include "atlas/output/Gmsh.h"
#include "atlas/util/CoordinateEnums.h"
#include "atlas/util/function/VortexRollup.h"

#include "pluto/pluto.h"

#include "tests/AtlasTestEnvironment.h"

using atlas::functionspace::NodeColumns;
using atlas::functionspace::StructuredColumns;
using atlas::util::Config;

namespace atlas {
namespace test {

//-----------------------------------------------------------------------------

static Config scheme() {
    Config scheme;
    scheme.set("type", "structured-cubic2D");
    scheme.set("halo", 2);
    scheme.set("name", "cubic");
    scheme.set("verbose",eckit::Resource<bool>("--verbose",false));
    scheme.set("sparse_matrix_multiply", "hicsparse");
    return scheme;
}

std::string input_gridname(const std::string& default_grid) {
    return eckit::Resource<std::string>("--input-grid", default_grid);
}

std::string output_gridname(const std::string& default_grid) {
    return eckit::Resource<std::string>("--output-grid", default_grid);
}

CASE("which scheme?") {
    Log::info() << scheme().getString("type") << std::endl;
}

template <typename Value>
struct AdjointTolerance {
    static const Value value;
};
template <>
const double AdjointTolerance<double>::value = 2.e-14;
template <>
const float AdjointTolerance<float>::value = 2.e-5;

template <typename Value, int rank>
void test_interpolation_structured_using_fs_API_for_fieldset_w_hicsparse_backend(const bool start_with_data_on_device) {
    auto pluto_trace = pluto::trace::enabled();
    pluto::trace::enable(false);

    Grid input_grid(input_gridname("O32"));
    Grid output_grid(output_gridname("O64"));

    // Cubic interpolation requires a StructuredColumns functionspace with 2 halos
    StructuredColumns input_fs(input_grid, scheme() | option::levels(rank == 1 ? 0 : 3));

    MeshGenerator meshgen("structured");
    Mesh output_mesh        = meshgen.generate(output_grid);
    FunctionSpace output_fs = NodeColumns{output_mesh, option::levels(rank == 1 ? 0 : 3)};

    auto lonlat = array::make_view<double, 2>(input_fs.xy());

    pluto::trace::enable(pluto_trace);
    FieldSet fields_source;
    FieldSet fields_target;
    for (idx_t f = 0; f < 3; ++f) {
        auto field_source = fields_source.add(input_fs.createField<Value>(option::name("source field " + std::to_string(f))));
        fields_target.add(output_fs.createField<Value>(option::name("target field " + std::to_string(f))));

        if (rank==1) {
            auto source = array::make_view<Value, 1>(field_source);
            for (idx_t n = 0; n < input_fs.size(); ++n) {
                source(n) = util::function::vortex_rollup(lonlat(n, LON), lonlat(n, LAT), 0.5);
            }
        }
        else if (rank == 2) {
            auto source = array::make_view<Value, 2>(field_source);
            for (idx_t n = 0; n < input_fs.size(); ++n) {
                for (idx_t k = 0; k < 3; ++k) {
                    source(n, k) = util::function::vortex_rollup(lonlat(n, LON), lonlat(n, LAT), 0.5 + double(k) / 2);
                }
            }
        }
        else {
            ATLAS_NOTIMPLEMENTED;
        }
    }

    if (start_with_data_on_device) {
        ATLAS_TRACE("Copy fields_source to device");
        fields_source.syncDevice();
    }

    ATLAS_TRACE_SCOPE("halo-exchange (fields_source)") {
        fields_source.haloExchange(start_with_data_on_device);
    }

    ATLAS_TRACE_SCOPE("with matrix") {
        Interpolation interpolation(scheme(), input_fs, output_fs);
        interpolation.execute(fields_source, fields_target);
        fields_target.syncHost();
        fields_target.haloExchange();
        ATLAS_TRACE_SCOPE("output") {
            output::Gmsh gmsh(scheme().getString("name") + "-multilevel-fieldset-output-with-matrix-" +
                              array::make_datatype<Value>().str() + ".msh",
                              Config("coordinates", "xy"));
            gmsh.write(output_mesh);
            gmsh.write(fields_target);
        }
    }

    ATLAS_TRACE_SCOPE("with matrix adjoint") {
        Interpolation interpolation(scheme() | Config("adjoint", true), input_fs, output_fs);

        std::vector<Value> AxAx(fields_source.field_names().size(), 0.);
        std::vector<Value> xAtAx(fields_source.field_names().size(), 0.);

        FieldSet fields_source_reference;
        for (const atlas::Field& field : fields_source) {
            Field field_ref(field.name(), field.datatype().kind(), field.shape());
            field_ref.set_levels(field.levels());

            if (rank == 1) {
                auto fieldInView  = array::make_view<const Value, 1>(field);
                auto fieldOutView = array::make_view<Value, 1>(field_ref);

                for (atlas::idx_t jn = 0; jn < field_ref.shape(0); ++jn) {
                        fieldOutView(jn) = fieldInView(jn);
                }
            }
            else if (rank == 2) {
                auto fieldInView  = array::make_view<const Value, 2>(field);
                auto fieldOutView = array::make_view<Value, 2>(field_ref);

                for (atlas::idx_t jn = 0; jn < field_ref.shape(0); ++jn) {
                    for (atlas::idx_t jl = 0; jl < field_ref.levels(); ++jl) {
                        fieldOutView(jn, jl) = fieldInView(jn, jl);
                    }
                }
            }
            else {
                ATLAS_NOTIMPLEMENTED;
            }
            fields_source_reference.add(field_ref);
        }

        interpolation.execute(fields_source, fields_target);
        fields_target.syncHost();

        std::size_t fIndx(0);
        auto source_names = fields_source.field_names();
        for (const std::string& s : fields_target.field_names()) {

            if (rank == 1) {
                auto target = array::make_view<const Value, 1>(fields_target[s]);
                auto source = array::make_view<const Value, 1>(fields_source[source_names[fIndx]]);

                for (idx_t n = 0; n < input_fs.size(); ++n) {
                    AxAx[fIndx] += source(n) * source(n);
                }

                for (idx_t n = 0; n < output_fs.size(); ++n) {
                    AxAx[fIndx] += target(n) * target(n);
                }
            }
            else if (rank == 2) {
                auto target = array::make_view<const Value, 2>(fields_target[s]);
                auto source = array::make_view<const Value, 2>(fields_source[source_names[fIndx]]);

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
            }
            else {
                ATLAS_NOTIMPLEMENTED;
            }

            fIndx += 1;
        }

        if (!start_with_data_on_device) {
            fields_source.syncHost();
            fields_target.syncHost();
            fields_source.deallocateDevice();
            fields_target.deallocateDevice();
        }

        interpolation.execute_adjoint(fields_source, fields_target);
        fields_source.syncHost();

        fIndx = 0;
        for (const std::string& s : fields_source.field_names()) {
            if (rank == 1) {
                auto source_reference = array::make_view<Value, 1>(fields_source_reference[s]);
                auto source           = array::make_view<Value, 1>(fields_source[s]);

                for (idx_t n = 0; n < input_fs.size(); ++n) {
                    xAtAx[fIndx] += source(n) * source_reference(n);
                }
            }
            else if (rank == 2) {
                auto source_reference = array::make_view<Value, 2>(fields_source_reference[s]);
                auto source           = array::make_view<Value, 2>(fields_source[s]);

                for (idx_t n = 0; n < input_fs.size(); ++n) {
                    for (idx_t k = 0; k < 3; ++k) {
                        xAtAx[fIndx] += source(n, k) * source_reference(n, k);
                    }
                }
            }
            else {
                ATLAS_NOTIMPLEMENTED;
            }
            fIndx += 1;
        }

        for (std::size_t t = 0; t < AxAx.size(); ++t) {
            Log::debug() << " Adjoint test t  = " << t << " (Ax).(Ax) = " << AxAx[t] << " x.(AtAx) = " << xAtAx[t]
                         << " std::abs( 1.0 - xAtAx[t]/AxAx[t] ) " << std::abs(1.0 - xAtAx[t] / AxAx[t])
                         << " AdjointTolerance<Value>::value " << AdjointTolerance<Value>::value << std::endl;

            EXPECT_APPROX_EQ(xAtAx[t] / AxAx[t] , 1.0, AdjointTolerance<Value>::value);
        }
    }
    pluto::trace::enable(pluto_trace);
}

constexpr bool start_on_host = false;
constexpr bool start_on_device = true;

CASE("test_interpolation_structured using fs API for fieldset with hicsparse backend (start with rank-2 double-precision data on host)") {
    test_interpolation_structured_using_fs_API_for_fieldset_w_hicsparse_backend<double,2>(start_on_host);
}

CASE("test_interpolation_structured using fs API for fieldset with hicsparse backend (start with rank-1 double-precision data on host)") {
    test_interpolation_structured_using_fs_API_for_fieldset_w_hicsparse_backend<double,1>(start_on_host);
}

CASE("test_interpolation_structured using fs API for fieldset with hicsparse backend (start with rank-2 double-precision data on device)") {
    test_interpolation_structured_using_fs_API_for_fieldset_w_hicsparse_backend<double,2>(start_on_device);
}

CASE("test_interpolation_structured using fs API for fieldset with hicsparse backend (start with rank-1 double-precision data on device)") {
    test_interpolation_structured_using_fs_API_for_fieldset_w_hicsparse_backend<double,1>(start_on_device);
}

#define SINGLE_PRECISION_WORKS 0
#if SINGLE_PRECISION_WORKS
CASE("test_interpolation_structured using fs API for fieldset with hicsparse backend (start with rank-2 single-precision data on host)") {
    test_interpolation_structured_using_fs_API_for_fieldset_w_hicsparse_backend<float,2>(start_on_host);
}

CASE("test_interpolation_structured using fs API for fieldset with hicsparse backend (start with rank-1 single-precision data on host)") {
    test_interpolation_structured_using_fs_API_for_fieldset_w_hicsparse_backend<float,1>(start_on_host);
}

CASE("test_interpolation_structured using fs API for fieldset with hicsparse backend (start with rank-2 single-precision data on device)") {
    test_interpolation_structured_using_fs_API_for_fieldset_w_hicsparse_backend<float,2>(start_on_device);
}

CASE("test_interpolation_structured using fs API for fieldset with hicsparse backend (start with rank-1 single-precision data on device)") {
    test_interpolation_structured_using_fs_API_for_fieldset_w_hicsparse_backend<float,1>(start_on_device);
}
#endif


}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
