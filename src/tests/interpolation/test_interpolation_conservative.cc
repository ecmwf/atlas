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
#include "atlas/util/function/SlottedCylinder.h"

#include "tests/AtlasTestEnvironment.h"


namespace atlas {
namespace test {

using ConservativeMethod = interpolation::method::ConservativeSphericalPolygonInterpolation;
using Statistics         = ConservativeMethod::Statistics;
using Metadata           = util::Metadata;
enum RemapStats {CONS = 0, CONS_MFREE, CONS2, CONS2_MFREE, CONS2_LIM, CONS2_MFREE_LIM, REMAPSTATS_SIZE};

void do_remapping_test(Grid src_grid, Grid tgt_grid, std::function<double(const PointLonLat&)> func,
                       std::vector<Metadata>& remap_stats, bool src_cell_data, bool tgt_cell_data) {
    std::string src_data_type = (src_cell_data ? "CellColumns(" : "NodeColumns(");
    std::string tgt_data_type = (tgt_cell_data ? "CellColumns(" : "NodeColumns(");
    Log::info() << "+-----------------------\n";
    Log::info() << src_data_type << src_grid.name() << ") --> " << tgt_data_type << tgt_grid.name() <<")\n";
    Log::info() << "+-----------------------\n";
    Log::info().indent();

    remap_stats.clear();
    remap_stats.resize(RemapStats::REMAPSTATS_SIZE);
    util::Config config("type", "conservative-spherical-polygon");
    config.set("order", 1);
    config.set("statistics.intersection", true);
    config.set("statistics.conservation", true);
    config.set("src_cell_data", src_cell_data);
    config.set("tgt_cell_data", tgt_cell_data);
    config.set("matrix_free", false);

    config.set("statistics.accuracy", true);
    auto conservative_interpolation = Interpolation(config, src_grid, tgt_grid);
    Log::info() << conservative_interpolation << std::endl << std::endl;

    const auto& src_fs = conservative_interpolation.source();
    const auto& tgt_fs = conservative_interpolation.target();
    auto src_field     = src_fs.createField<double>();
    auto tgt_field     = tgt_fs.createField<double>();

    ConservativeMethod& consMethod = dynamic_cast<ConservativeMethod&>(*conservative_interpolation.get());
    ATLAS_TRACE_SCOPE("initial condition") {
        auto src_vals = array::make_view<double, 1>(src_field);
        for (idx_t spt = 0; spt < src_vals.size(); ++spt) {
            auto p = consMethod.src_points(spt);
            PointLonLat pll;
            eckit::geometry::Sphere::convertCartesianToSpherical(1., p, pll);
            src_vals(spt) = func(pll);
        }
    }

    ATLAS_TRACE_SCOPE("1st order projection matrix-version") {
        remap_stats[RemapStats::CONS] = conservative_interpolation.execute(src_field, tgt_field);
        tgt_field.haloExchange();
        consMethod.statistics().compute_accuracy(conservative_interpolation, tgt_field, func, &remap_stats[RemapStats::CONS]);
    }
    ATLAS_TRACE_SCOPE("1st order projection matrix-free version") {
        config.set("matrix_free", true);
        config.set("statistics.accuracy", true);
        auto conservative_interpolation = Interpolation(config, src_grid, tgt_grid);
        remap_stats[RemapStats::CONS_MFREE] = conservative_interpolation.execute(src_field, tgt_field);
        tgt_field.haloExchange();
        consMethod.statistics().compute_accuracy(conservative_interpolation, tgt_field, func, &remap_stats[RemapStats::CONS_MFREE]);
    }
    ATLAS_TRACE_SCOPE("2nd order projection matrix-version") {
        config.set("order", 2);
        config.set("matrix_free", false);
        config.set("statistics.accuracy", true);
        if (src_cell_data && tgt_cell_data) {
            config.set("limiter", "none");
        }
        conservative_interpolation = Interpolation(config, src_grid, tgt_grid);
        Log::info() << conservative_interpolation << std::endl;
        remap_stats[RemapStats::CONS2] = conservative_interpolation.execute(src_field, tgt_field);
        auto& consMethod_2 = dynamic_cast<ConservativeMethod&>(*conservative_interpolation.get());
        tgt_field.haloExchange();
        consMethod_2.statistics().compute_accuracy(conservative_interpolation, tgt_field, func, &remap_stats[RemapStats::CONS2]);
        if (src_cell_data && tgt_cell_data) {
            config.set("limiter", "zeroslope");
            conservative_interpolation = Interpolation(config, src_grid, tgt_grid);
            Log::info() << conservative_interpolation << std::endl;
            remap_stats[RemapStats::CONS2_LIM] = conservative_interpolation.execute(src_field, tgt_field);
            auto& consMethod_2 = dynamic_cast<ConservativeMethod&>(*conservative_interpolation.get());
            tgt_field.haloExchange();
            consMethod_2.statistics().compute_accuracy(conservative_interpolation, tgt_field, func, &remap_stats[RemapStats::CONS2_LIM]);
        }
        else {
            // remap_stat[RemapStats::CONS2_LIM].reset
        }
    }
    ATLAS_TRACE_SCOPE("2nd order projection matrix-free version") {
        config.set("order", 2);
        config.set("matrix_free", true);
        config.set("statistics.accuracy", true);
        conservative_interpolation = Interpolation(config, src_grid, tgt_grid);
        Log::info() << conservative_interpolation << std::endl;
        remap_stats[RemapStats::CONS2_MFREE] = conservative_interpolation.execute(src_field, tgt_field);
        auto& consMethod_2 = dynamic_cast<ConservativeMethod&>(*conservative_interpolation.get());
        tgt_field.haloExchange();
        consMethod_2.statistics().compute_accuracy(conservative_interpolation, tgt_field, func, &remap_stats[RemapStats::CONS2_MFREE]);
        if (src_cell_data && tgt_cell_data) {
            config.set("limiter", "zeroslope");
            conservative_interpolation = Interpolation(config, src_grid, tgt_grid);
            Log::info() << conservative_interpolation << std::endl;
            remap_stats[RemapStats::CONS2_MFREE_LIM] = conservative_interpolation.execute(src_field, tgt_field);
            auto& consMethod_2 = dynamic_cast<ConservativeMethod&>(*conservative_interpolation.get());
            tgt_field.haloExchange();
            consMethod_2.statistics().compute_accuracy(conservative_interpolation, tgt_field, func, &remap_stats[RemapStats::CONS2_MFREE_LIM]);
        }
        else {
            remap_stats[RemapStats::CONS2_MFREE_LIM].set("errors.to_exact_solution_sum", -1.);
            remap_stats[RemapStats::CONS2_MFREE_LIM].set("errors.to_exact_solution_max", -1.);
        }
    }

    ATLAS_TRACE_SCOPE("test caching") {
        // We can create the interpolation without polygon intersections
        auto cache = interpolation::Cache(conservative_interpolation);
        // cache = ConservativeMethod::Cache + MatrixCache (1st order)
        util::Config cfg(option::type("conservative-spherical-polygon"));
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
        {
            ATLAS_TRACE("cached -> 2nd order matrix-free");
            cfg.set("matrix_free", true);
            cfg.set("order", 2);
            auto interpolation = Interpolation(cfg, src_grid, tgt_grid, cache);
            Log::info() << interpolation << std::endl;
            interpolation.execute(src_field, tgt_field);
            Log::info() << std::endl;
        }
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
}

void check(const std::vector<Metadata>& remap_stat, std::array<double, 6> tol, bool cell_data) {
    double err_mv;
    double err_mf;
    double err_mv_lim;
    double err_mf_lim;

    remap_stat[RemapStats::CONS].get("errors.sum_src_areas_minus_sum_tgt_areas", err_mv);
    Log::info() << "Polygon area computation (new < ref)                 =  (" << err_mv << " < " << tol[0] << ")" << std::endl;
    EXPECT(err_mv < tol[0]);
    remap_stat[RemapStats::CONS].get("errors.intersections_covering_tgt_cells_sum", err_mv);
    Log::info() << "Polygon intersection (new < ref)                     =  (" << err_mv << " < " << tol[1] << ")" << std::endl;
    EXPECT(err_mv < tol[1]);

    remap_stat[RemapStats::CONS].get("errors.to_exact_solution_sum", err_mv);
    Log::info() << "\n1st order accuracy (new < ref)                       =  (" << std::abs(err_mv) << " < " << tol[2] << ")" << std::endl;
    EXPECT(std::abs(err_mv) < tol[2]);
    remap_stat[RemapStats::CONS_MFREE].get("errors.to_exact_solution_sum", err_mf);
    Log::info() << "1st order matrix-free accuracy (new < ref)           =  (" << std::abs(err_mf) << " < " << tol[2] << ")" << std::endl;
    Log::info() << "                    |matrix - matrix_free|           =   " << std::abs(err_mv - err_mf) << std::endl;
    EXPECT(std::abs(err_mf) < tol[2]);

    remap_stat[RemapStats::CONS2].get("errors.to_exact_solution_sum", err_mv);
    Log::info() << "\n2nd order accuracy (new < ref)                       =  (" << std::abs(err_mv) << " < " << tol[3] << ")" << std::endl;
    EXPECT(std::abs(err_mv) < tol[3]);
    remap_stat[RemapStats::CONS2_MFREE].get("errors.to_exact_solution_sum", err_mf);
    Log::info() << "2nd order matrix-free accuracy (new < ref)           =  (" << std::abs(err_mf) << " < " << tol[3] << ")" << std::endl;
    Log::info() << "                    |matrix - matrix_free|           =   " << std::abs(err_mv - err_mf) << std::endl;
    EXPECT(std::abs(err_mf) < tol[3]);
    if (cell_data) {
        remap_stat[RemapStats::CONS2_LIM].get("errors.to_exact_solution_sum", err_mv_lim);
        Log::info() << "2nd order lim. matrix-ver. accuracy (new < ref)      =  (" << std::abs(err_mv_lim) << " < " << tol[3] << ")" << std::endl;
        Log::info() << "          |lim. matrix-ver. - matrix_ver.|           =   " << std::abs(err_mv_lim - err_mv) << std::endl;
        remap_stat[RemapStats::CONS2_MFREE_LIM].get("errors.to_exact_solution_sum", err_mf_lim);
        Log::info() << "2nd order lim. matrix-free accuracy (new < ref)      =  (" << std::abs(err_mf_lim) << " < " << tol[3] << ")" << std::endl;
        Log::info() << "          |lim. matrix_free - matrix_free|           =   " << std::abs(err_mf - err_mf_lim) << std::endl;
    }

    remap_stat[RemapStats::CONS].get("errors.conservation", err_mv);
    Log::info() << "\n1st order conservation (new < ref)                   =  (" << std::abs(err_mv) << " < " << tol[4] << ")" << std::endl;
    EXPECT(std::abs(err_mv) < tol[4]);
    remap_stat[RemapStats::CONS_MFREE].get("errors.conservation", err_mf);
    Log::info() << "1st order matrix-free conservation (new < ref)       =  (" << std::abs(err_mf) << " < " << tol[4] << ")" << std::endl;
    EXPECT(std::abs(err_mf) < tol[4]);

    remap_stat[RemapStats::CONS2].get("errors.conservation", err_mv);
    Log::info() << "\n2nd order conservation (new < ref)                   =  (" << std::abs(err_mv) << " < " << tol[5] << ")" << std::endl;
    EXPECT(std::abs(err_mv) < tol[5]);remap_stat[RemapStats::CONS2_MFREE].get("errors.conservation", err_mf);
    Log::info() << "2nd order matrix-free conservation (new < ref)       =  (" << std::abs(err_mf) << " < " << tol[5] << ")" << std::endl;
    EXPECT(std::abs(err_mf) < tol[5]);
    if (cell_data) {
        remap_stat[RemapStats::CONS2_LIM].get("errors.conservation", err_mv_lim);
        Log::info() << "2nd order lim. matrix-ver. conservation (new < ref)  =  (" << std::abs(err_mv_lim) << " < " << tol[5] << ")" << std::endl;
        remap_stat[RemapStats::CONS2_MFREE_LIM].get("errors.conservation", err_mf_lim);
        Log::info() << "2nd order lim. matrix-free conservation (new < ref)  =  (" << std::abs(err_mf_lim) << " < " << tol[5] << ")" << std::endl;
    }
    Log::info().unindent();
}

CASE("test_interpolation_conservative") {
    std::vector<Metadata> remap_stats(RemapStats::REMAPSTATS_SIZE);

    SECTION("analytic constfunc") {
        auto func = [](const PointLonLat& p) { return 1.; };
        bool src_cell_data = true;
        bool tgt_cell_data = true;
        do_remapping_test(Grid("O32"), Grid("H12"), func, remap_stats, src_cell_data, tgt_cell_data);
        check(remap_stats, {1.0e-13, 1.0e-13, 1.0e-13, 1.0e-13, 1.0e-13, 1.0e-13}, src_cell_data && tgt_cell_data);
    }

    SECTION("vortex_rollup") {
        auto func = [](const PointLonLat& p) {
            return util::function::vortex_rollup(p[0], p[1], 0.5);
        };

        bool src_cell_data = true;
        bool tgt_cell_data = true;
        do_remapping_test(Grid("O16"), Grid("H12"), func, remap_stats, src_cell_data, tgt_cell_data);
        check(remap_stats, {2.0e-14, 1.0e-14, 0.0051927, 0.0025274, 1.0e-15, 1.5e-08}, src_cell_data && tgt_cell_data);

        src_cell_data = true;
        tgt_cell_data = false;
        do_remapping_test(Grid("O16"), Grid("H12"), func, remap_stats, src_cell_data, tgt_cell_data);
        check(remap_stats, {2.0e-14, 1.0e-14, 0.0054418, 0.0028356, 1.0e-15, 3.0e-09}, src_cell_data && tgt_cell_data);

        src_cell_data = false;
        tgt_cell_data = true;
        do_remapping_test(Grid("O16"), Grid("H12"), func, remap_stats, src_cell_data, tgt_cell_data);
        check(remap_stats, {2.0e-14, 1.0e-14, 0.0062701, 0.0029492, 5.0e-16, 6.0e-10}, src_cell_data && tgt_cell_data);

        src_cell_data = false;
        tgt_cell_data = false;
        do_remapping_test(Grid("O16"), Grid("H12"), func, remap_stats, src_cell_data, tgt_cell_data);
        check(remap_stats, {2.0e-14, 1.0e-14, 0.0064164, 0.0030295, 5.0e-16, 5.0e-13}, src_cell_data && tgt_cell_data);
    }

    SECTION("slotted_cylinder") {
        auto func = [](const PointLonLat& p) {
            return util::function::SlottedCylinder(p[0], p[1]);
        };

        bool src_cell_data = true;
        bool tgt_cell_data = true;
        do_remapping_test(Grid("O16"), Grid("H12"), func, remap_stats, src_cell_data, tgt_cell_data);
        check(remap_stats, {2.0e-14, 1.0e-14, 0.0624738, 0.0620427, 1.0e-15, 5.0e-08}, src_cell_data && tgt_cell_data);
    }
}

}  // namespace test
}  // namespace atlas


int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
