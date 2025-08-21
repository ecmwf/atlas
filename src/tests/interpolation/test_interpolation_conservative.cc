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

void do_remapping_test(Grid src_grid, Grid tgt_grid, std::function<double(const PointLonLat&)> func,
                       Statistics& remap_stat_1, Statistics& remap_stat_2, bool src_cell_data, bool tgt_cell_data) {
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

    {
        ATLAS_TRACE("initial condition");
        // A bit of a hack here...
        ConservativeMethod& consMethod = dynamic_cast<ConservativeMethod&>(*conservative_interpolation.get());
        for (idx_t spt = 0; spt < src_vals.size(); ++spt) {
            auto p = consMethod.src_points(spt);
            PointLonLat pll;
            eckit::geometry::Sphere::convertCartesianToSpherical(1., p, pll);
            src_vals(spt) = func(pll);
        }
    }

    // project source field to target mesh in 1st order
    remap_stat_1 = conservative_interpolation.execute(src_field, tgt_field);
    remap_stat_1.accuracy(conservative_interpolation, tgt_field, func);

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
        if (src_cell_data and tgt_cell_data) {
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
        if (src_cell_data and tgt_cell_data) {
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

    {
        // project source field to target mesh in 2nd order
        config.set("order", 2);
        conservative_interpolation = Interpolation(config, src_grid, tgt_grid);
        Log::info() << conservative_interpolation << std::endl;
        remap_stat_2 = conservative_interpolation.execute(src_field, tgt_field);
        remap_stat_2.accuracy(conservative_interpolation, tgt_field, func);
    }
}

void check(const Statistics remap_stat_1, Statistics remap_stat_2, std::array<double, 6> tol) {
    double err;
    // check polygon intersections
    err = remap_stat_1.errors[Statistics::Errors::SRCTGT_INTERSECTPLG_DIFF];
    Log::info() << "Polygon area computation (new < ref) =  (" << err << " < " << tol[0] << ")" << std::endl;
    EXPECT(err < tol[0]);
    err = remap_stat_1.errors[Statistics::Errors::TGT_INTERSECTPLG_L1];
    Log::info() << "Polygon intersection (new < ref) =  (" << err << " < " << tol[1] << ")" << std::endl;
    EXPECT(err < tol[1]);

    // check remap accuracy
    err = std::abs(remap_stat_1.errors[Statistics::Errors::REMAP_L2]);
    Log::info() << "1st order accuracy (new < ref) =  (" << err << " < " << tol[2] << ")" << std::endl;
    EXPECT(err < tol[2]);
    err = std::abs(remap_stat_2.errors[Statistics::Errors::REMAP_L2]);
    Log::info() << "2nd order accuracy (new < ref) =  (" << err << " < " << tol[3] << ")" << std::endl;
    EXPECT(err < tol[3]);

    // check mass conservation
    err = std::abs(remap_stat_1.errors[Statistics::Errors::REMAP_CONS]);
    Log::info() << "1st order conservation (new < ref) =  (" << err << " < " << tol[4] << ")" << std::endl;
    EXPECT(err < tol[4]);
    err = std::abs(remap_stat_2.errors[Statistics::Errors::REMAP_CONS]);
    Log::info() << "2nd order conservation (new < ref) =  (" << err << " < " << tol[5] << ")" << std::endl;
    EXPECT(err < tol[5]);
    Log::info().unindent();
}

CASE("test_interpolation_conservative") {
    SECTION("analytic constfunc") {
        auto func = [](const PointLonLat& p) { return 1.; };
        Statistics remap_stat_1;
        Statistics remap_stat_2;
        bool src_cell_data = true;
        bool tgt_cell_data = true;
        do_remapping_test(Grid("O32"), Grid("H12"), func, remap_stat_1, remap_stat_2, src_cell_data, tgt_cell_data);
        check(remap_stat_1, remap_stat_2, {1.0e-13, 1.0e-13, 1.0e-13, 1.0e-13, 1.0e-13, 1.0e-13});
    }

    SECTION("vortex_rollup") {
        auto func = [](const PointLonLat& p) {
            return util::function::vortex_rollup(p[0], p[1], 0.5);
        };
        Statistics remap_stat_1;
        Statistics remap_stat_2;

        bool src_cell_data = true;
        bool tgt_cell_data = true;
        do_remapping_test(Grid("O16"), Grid("H12"), func, remap_stat_1, remap_stat_2, src_cell_data, tgt_cell_data);
        check(remap_stat_1, remap_stat_2, {1.0e-13, 1.0e-12, 7.0e-3, 4.0e-3, 1.0e-15, 1.5e-8});

        src_cell_data = true;
        tgt_cell_data = false;
        do_remapping_test(Grid("O16"), Grid("H12"), func, remap_stat_1, remap_stat_2, src_cell_data, tgt_cell_data);
        check(remap_stat_1, remap_stat_2, {1.0e-13, 1.0e-12, 7.0e-3, 4.0e-3, 1.0e-15, 1.0e-8});

        src_cell_data = false;
        tgt_cell_data = true;
        do_remapping_test(Grid("O16"), Grid("H12"), func, remap_stat_1, remap_stat_2, src_cell_data, tgt_cell_data);
        check(remap_stat_1, remap_stat_2, {1.0e-13, 1.0e-12, 7.0e-3, 4.0e-3, 1.0e-15, 2.0e-8});

        src_cell_data = false;
        tgt_cell_data = false;
        do_remapping_test(Grid("O16"), Grid("H12"), func, remap_stat_1, remap_stat_2, src_cell_data, tgt_cell_data);
        check(remap_stat_1, remap_stat_2, {1.0e-12, 1.0e-12, 7.0e-3, 4.0e-3, 1.0e-15, 1.0e-8});
    }
}

}  // namespace test
}  // namespace atlas


int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
