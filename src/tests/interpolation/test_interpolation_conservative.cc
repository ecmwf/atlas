/*
 * (C) Copyright 1996- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */


#include <cmath>

#include "eckit/geometry/Sphere.h"
#include "eckit/types/FloatCompare.h"

#include "atlas/array.h"
#include "atlas/array/MakeView.h"
#include "atlas/field.h"
#include "atlas/grid.h"
#include "atlas/interpolation.h"
#include "atlas/interpolation/method/unstructured/ConservativeSphericalPolygonInterpolation.h"
#include "atlas/mesh.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/meshgenerator.h"
#include "atlas/option.h"
#include "atlas/util/Config.h"
#include "atlas/util/function/VortexRollup.h"

#include "tests/AtlasTestEnvironment.h"


namespace atlas {
namespace test {

using ConservativeMethod = interpolation::method::ConservativeSphericalPolygonInterpolation;
using Statistics         = ConservativeMethod::Statistics;
using Metadata           = util::Metadata;

void do_remapping_test(Grid src_grid, Grid tgt_grid, std::function<double(const PointLonLat&)> func,
                       Metadata& remap_stat_1, Metadata& remap_stat_2, bool src_cell_data, bool tgt_cell_data) {
    std::string src_data_type = (src_cell_data ? "CellColumns(" : "NodeColumns(");
    std::string tgt_data_type = (tgt_cell_data ? "CellColumns(" : "NodeColumns(");
    Log::info() << "+-----------------------\n";
    Log::info() << src_data_type << src_grid.name() << ") --> " << tgt_data_type << tgt_grid.name() <<")\n";
    Log::info() << "+-----------------------\n";
    Log::info().indent();

    // setup conservative remap: compute weights, polygon intersection, etc
    util::Config config("type", "conservative-spherical-polygon");
    config.set("order", 1);
    config.set("validate", true);
    config.set("statistics.intersection", true);
    config.set("statistics.conservation", true);
    config.set("src_cell_data", src_cell_data);
    config.set("tgt_cell_data", tgt_cell_data);

    auto conservative_interpolation = Interpolation(config, src_grid, tgt_grid);
    Log::info() << conservative_interpolation << std::endl;
    Log::info() << std::endl;

    // create source field from analytic function "func"
    const auto& src_fs = conservative_interpolation.source();
    const auto& tgt_fs = conservative_interpolation.target();
    auto src_field     = src_fs.createField<double>();
    auto tgt_field     = tgt_fs.createField<double>();
    auto src_vals      = array::make_view<double, 1>(src_field);

    // A bit of a hack here...
    ConservativeMethod& consMethod = dynamic_cast<ConservativeMethod&>(*conservative_interpolation.get());
    {
        ATLAS_TRACE("initial condition");

        for (idx_t spt = 0; spt < src_vals.size(); ++spt) {
            auto p = consMethod.src_points(spt);
            PointLonLat pll;
            eckit::geometry::Sphere::convertCartesianToSpherical(1., p, pll);
            src_vals(spt) = func(pll);
        }
    }

    remap_stat_1 = conservative_interpolation.execute(src_field, tgt_field);
    tgt_field.haloExchange();
    consMethod.statistics().compute_accuracy(conservative_interpolation, tgt_field, func, &remap_stat_1);

    ATLAS_TRACE_SCOPE("test caching") {
        // We can create the interpolation without polygon intersections
        auto cache = interpolation::Cache(conservative_interpolation);
        // cache = ConservativeMethod::Cache + MatrixCache (1st order)
        util::Config cfg(option::type("conservative-spherical-polygon"));
        config.set("validate", true);
        {
            ATLAS_TRACE("cached -> 1st order using cached matrix");
            cfg.set("matrix_free", false);
            cfg.set("order", 1);
            auto interpolation = Interpolation(cfg, src_grid, tgt_grid, cache);
            Log::info() << interpolation << std::endl;
            interpolation.execute(src_field, tgt_field);
            Log::info() << std::endl;
        }
        {
            ATLAS_TRACE("cached -> 1st order constructing new matrix");
            cfg.set("matrix_free", false);
            cfg.set("order", 1);
            auto cache_without_matrix =
                ConservativeMethod::Cache(cache);  // to mimick when cache was created with matrix_free option
            auto interpolation = Interpolation(cfg, src_grid, tgt_grid, cache_without_matrix);
            Log::info() << interpolation << std::endl;
            interpolation.execute(src_field, tgt_field);
            Log::info() << std::endl;
        }
        {
            ATLAS_TRACE("cached -> 1st order matrix-free");
            cfg.set("matrix_free", true);
            cfg.set("order", 1);
            auto interpolation = Interpolation(cfg, src_grid, tgt_grid, cache);
            Log::info() << interpolation << std::endl;
            interpolation.execute(src_field, tgt_field);
            Log::info() << std::endl;
        }
        auto cache_2 = interpolation::Cache{};
        {
            ATLAS_TRACE("cached -> 2nd order constructing new matrix");
            cfg.set("matrix_free", false);
            cfg.set("order", 2);
            auto interpolation = Interpolation(cfg, src_grid, tgt_grid, cache);
            Log::info() << interpolation << std::endl;
            interpolation.execute(src_field, tgt_field);
            cache_2 = interpolation.createCache();
            Log::info() << std::endl;
        }
        // TODO: With the PR 318 we switch to interating over target elements, which in return requires a lot of changes.
        // The following code requires reimplementation of the matrix-free 2nd order method which will come after this PR.
        // Hence, we temporary disable this code.
        //
        // if (src_cell_data and tgt_cell_data) {
        //     ATLAS_TRACE("cached -> 2nd order matrix-free");
        //     cfg.set("matrix_free", true);
        //     cfg.set("order", 2);
        //     auto interpolation = Interpolation(cfg, src_grid, tgt_grid, cache);
        //     Log::info() << interpolation << std::endl;
        //     interpolation.execute(src_field, tgt_field);
        //     Log::info() << std::endl;
        // }
        {
            ATLAS_TRACE("cached -> 2nd order using cached matrix");
            cfg.set("matrix_free", false);
            cfg.set("order", 2);
            auto interpolation = Interpolation(cfg, src_grid, tgt_grid, cache_2);
            Log::info() << interpolation << std::endl;
            interpolation.execute(src_field, tgt_field);
            Log::info() << std::endl;
        }
    }

    {
        ATLAS_TRACE("2nd order projection");
        config.set("order", 2);
        conservative_interpolation = Interpolation(config, src_grid, tgt_grid);
        Log::info() << conservative_interpolation << std::endl;
        remap_stat_2 = conservative_interpolation.execute(src_field, tgt_field);
        auto& consMethod_2 = dynamic_cast<ConservativeMethod&>(*conservative_interpolation.get());
        tgt_field.haloExchange();
        consMethod_2.statistics().compute_accuracy(conservative_interpolation, tgt_field, func, &remap_stat_2);
    }
}

void check(const Metadata remap_stat_1, Metadata remap_stat_2, std::array<double, 6> tol) {
    double err;
    // check polygon intersections
    remap_stat_1.get("errors.sum_src_areas_minus_sum_tgt_areas", err);
    Log::info() << "Polygon area computation (new < ref) =  (" << err << " < " << tol[0] << ")" << std::endl;
    EXPECT(err < tol[0]);
    remap_stat_1.get("errors.intersections_covering_tgt_cells_sum", err);
    Log::info() << "Polygon intersection (new < ref) =  (" << err << " < " << tol[1] << ")" << std::endl;
    EXPECT(err < tol[1]);

    // check remap accuracy
    remap_stat_1.get("errors.to_solution_sum", err);
    Log::info() << "1st order accuracy (new < ref) =  (" << std::abs(err) << " < " << tol[2] << ")" << std::endl;
    EXPECT(std::abs(err) < tol[2]);
    remap_stat_2.get("errors.to_solution_sum", err);
    Log::info() << "2nd order accuracy (new < ref) =  (" << std::abs(err) << " < " << tol[3] << ")" << std::endl;
    EXPECT(std::abs(err) < tol[3]);

    // check mass conservation
    remap_stat_1.get("errors.conservation_error", err);
    Log::info() << "1st order conservation (new < ref) =  (" << std::abs(err) << " < " << tol[4] << ")" << std::endl;
    EXPECT(std::abs(err) < tol[4]);
    remap_stat_2.get("errors.conservation_error", err);
    Log::info() << "2nd order conservation (new < ref) =  (" << std::abs(err) << " < " << tol[5] << ")" << std::endl;
    EXPECT(std::abs(err) < tol[5]);
    Log::info().unindent();
}

CASE("test_interpolation_conservative") {
    SECTION("analytic constfunc") {
        auto func = [](const PointLonLat& p) { return 1.; };
        Metadata remap_stat_1;
        Metadata remap_stat_2;
        bool src_cell_data = true;
        bool tgt_cell_data = true;
        do_remapping_test(Grid("O32"), Grid("H12"), func, remap_stat_1, remap_stat_2, src_cell_data, tgt_cell_data);
        check(remap_stat_1, remap_stat_2, {1.0e-13, 1.0e-13, 1.0e-13, 1.0e-13, 1.0e-13, 1.0e-13});
    }

    SECTION("vortex_rollup") {
        auto func = [](const PointLonLat& p) {
            return util::function::vortex_rollup(p[0], p[1], 0.5);
        };
        Metadata remap_stat_1;
        Metadata remap_stat_2;

        bool src_cell_data = true;
        bool tgt_cell_data = true;
        do_remapping_test(Grid("O16"), Grid("H12"), func, remap_stat_1, remap_stat_2, src_cell_data, tgt_cell_data);
        check(remap_stat_1, remap_stat_2, {1.0e-13, 1.0e-12, 0.0051927, 0.0025275, 1.0e-15, 1.5e-08});

        src_cell_data = true;
        tgt_cell_data = false;
        do_remapping_test(Grid("O16"), Grid("H12"), func, remap_stat_1, remap_stat_2, src_cell_data, tgt_cell_data);
        check(remap_stat_1, remap_stat_2, {1.0e-13, 1.0e-12, 0.0054418, 0.0028355, 1.0e-15, 5.0e-09});

        src_cell_data = false;
        tgt_cell_data = true;
        do_remapping_test(Grid("O16"), Grid("H12"), func, remap_stat_1, remap_stat_2, src_cell_data, tgt_cell_data);
        check(remap_stat_1, remap_stat_2, {1.0e-13, 1.0e-12, 0.0062715, 0.0029492, 1.0e-15, 2.0e-09});

        src_cell_data = false;
        tgt_cell_data = false;
        do_remapping_test(Grid("O16"), Grid("H12"), func, remap_stat_1, remap_stat_2, src_cell_data, tgt_cell_data);
        check(remap_stat_1, remap_stat_2, {1.0e-12, 1.0e-12, 0.0064164, 0.0030295, 1.0e-15, 1.0e-12});
    }
}

}  // namespace test
}  // namespace atlas


int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
