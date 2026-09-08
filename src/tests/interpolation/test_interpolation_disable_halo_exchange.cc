/*
 * (C) Crown Copyright 2026 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */
#include <eckit/testing/Test.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>

#include "atlas/functionspace.h"
#include "atlas/functionspace/StructuredColumns.h"
#include "atlas/grid.h"
#include "atlas/grid/Partitioner.h"
#include "atlas/interpolation.h"

#include "atlas/interpolation/Interpolation.h"
#include "atlas/util/Config.h"
#include "tests/AtlasTestEnvironment.h"

namespace atlas {
namespace test {

//-----------------------------------------------------------------------------

class FunctionSpaces {
public:
    FunctionSpaces():
        source_functionspace_{Grid("O8"), option::halo(1)},
        target_functionspace_{Grid("O24"), grid::MatchingPartitioner(source_functionspace_)} {}
    const functionspace::StructuredColumns& source() const { return source_functionspace_; }
    const functionspace::StructuredColumns& target() const { return target_functionspace_; }

private:
    functionspace::StructuredColumns source_functionspace_{};
    functionspace::StructuredColumns target_functionspace_{};
};

std::pair<Field, Field> make_fields(const FunctionSpaces& fs) {
    const auto source_field = fs.source().createField<double>(option::name("source"));
    const auto target_field = fs.target().createField<double>(option::name("target"));
    return std::make_pair(source_field, target_field);
}

void assign_field(Field field, double owned_value, double ghost_value) {
    const auto ghost = array::make_view<int, 1>(field.functionspace().ghost());
    auto view        = array::make_view<double, 1>(field);
    for (int i = 0; i < view.size(); ++i) {
        view(i) = (ghost(i) ? ghost_value : owned_value);
    }
}

Interpolation create_interp_object(const FunctionSpaces& fs, bool perform_halo_exchange = true) {
    auto conf = option::type("structured-bilinear") | util::Config("ajoint", true);

    if (!perform_halo_exchange) {
        conf.set("perform_halo_exchange", false);
    }

    return Interpolation(conf, fs.source(), fs.target());
}


CASE("Forward interpolation") {
    if constexpr (std::numeric_limits<double>::has_quiet_NaN) {
        const auto function_spaces = FunctionSpaces();

        const auto make_and_init_fields = [&] {
            auto [source_field, target_field] = make_fields(function_spaces);
            assign_field(source_field, 1., std::numeric_limits<double>::quiet_NaN());
            assign_field(target_field, 0., 0.);
            return std::make_pair(source_field, target_field);
        };

        const auto has_nan = [](const Field& field) {
            const auto view = array::make_view<const double, 1>(field);
            return std::any_of(view.data(), view.data() + view.size(), [](double value) { return std::isnan(value); });
        };

        SECTION("Automatic halo exchange enabled (default behaviour)") {
            auto [source_field, target_field] = make_and_init_fields();
            const auto interp                 = create_interp_object(function_spaces);
            interp.execute(source_field, target_field);
            EXPECT(!source_field.dirty());
            EXPECT(!has_nan(target_field));
        }
        SECTION("Automatic halo exchange disabled") {
            auto [source_field, target_field] = make_and_init_fields();
            const auto interp                 = create_interp_object(function_spaces, false);
            interp.execute(source_field, target_field);
            EXPECT(source_field.dirty());
            EXPECT(has_nan(target_field));
        }
        SECTION("Manual halo exchange") {
            auto [source_field, target_field] = make_and_init_fields();
            source_field.haloExchange();
            const auto interp = create_interp_object(function_spaces, false);
            interp.execute(source_field, target_field);
            EXPECT(!source_field.dirty());
            EXPECT(!has_nan(target_field));
        }
    }
}

CASE("Adjoint interpolation") {
    const auto function_spaces = FunctionSpaces();

    const auto make_and_init_fields = [&] {
        auto [source_field, target_field] = make_fields(function_spaces);
        assign_field(source_field, 0., 0.);
        assign_field(target_field, 1., 1.);
        return std::make_pair(source_field, target_field);
    };

    const auto has_empty_halos = [](Field& field) {
        const auto ghost = array::make_view<int, 1>(field.functionspace().ghost());
        const auto view  = array::make_view<double, 1>(field);
        for (int i = 0; i < view.size(); ++i) {
            if (ghost(i) && view(i) != 0.) {
                return false;
            }
        }
        return true;
    };

    SECTION("Automatic halo exchange enabled (default behaviour)") {
        auto [source_field, target_field] = make_and_init_fields();
        const auto interp                 = create_interp_object(function_spaces);
        interp.execute_adjoint(source_field, target_field);
        EXPECT(!source_field.dirty());
        EXPECT(has_empty_halos(source_field));
    }
    SECTION("Automatic halo exchange disabled") {
        auto [source_field, target_field] = make_and_init_fields();
        const auto interp                 = create_interp_object(function_spaces, false);
        interp.execute_adjoint(source_field, target_field);
        EXPECT(source_field.dirty());
        EXPECT(!has_empty_halos(source_field));
    }
    SECTION("Manual halo exchange") {
        auto [source_field, target_field] = make_and_init_fields();
        const auto interp                 = create_interp_object(function_spaces, false);
        interp.execute_adjoint(source_field, target_field);
        source_field.adjointHaloExchange();
        EXPECT(!source_field.dirty());
        EXPECT(has_empty_halos(source_field));
    }
}
//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
