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

#include "tests/AtlasTestEnvironment.h"


namespace atlas {
namespace test {

using ConservativeMethod = interpolation::method::ConservativeSphericalPolygonInterpolation;
using Statistics         = ConservativeMethod::Statistics;

void do_remapping_test(Grid src_grid, Grid tgt_grid, std::function<double(const PointLonLat&)> func,
                       Statistics& remap_stat_1, Statistics& remap_stat_2) {
    Log::info().indent();
    // setup conservative remap: compute weights, polygon intersection, etc
    util::Config config("type", "conservative-spherical-polygon");
    config.set("order", 1);
    config.set("validate", true);
    config.set("statistics.intersection", true);
    config.set("statistics.conservation", true);

    auto conservative_interpolation = Interpolation(config, src_grid, tgt_grid);
    Log::info() << conservative_interpolation << std::endl;

    // create source field from analytic function "func"
    const auto& src_fs = conservative_interpolation.source();
    const auto& tgt_fs = conservative_interpolation.target();
    auto src_field     = src_fs.createField<double>();
    auto tgt_field     = tgt_fs.createField<double>();
    auto src_vals      = array::make_view<double, 1>(src_field);
    auto tgt_vals      = array::make_view<double, 1>(tgt_field);

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
        }
        {
            ATLAS_TRACE("cached -> 1st order matrix-free");
            cfg.set("matrix_free", true);
            cfg.set("order", 1);
            auto interpolation = Interpolation(cfg, src_grid, tgt_grid, cache);
            Log::info() << interpolation << std::endl;
            interpolation.execute(src_field, tgt_field);
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
        }
        {
            ATLAS_TRACE("cached -> 2nd order matrix-free");
            cfg.set("matrix_free", true);
            cfg.set("order", 2);
            auto interpolation = Interpolation(cfg, src_grid, tgt_grid, cache);
            Log::info() << interpolation << std::endl;
            interpolation.execute(src_field, tgt_field);
        }
        {
            ATLAS_TRACE("cached -> 2nd order using cached matrix");
            cfg.set("matrix_free", false);
            cfg.set("order", 2);
            auto interpolation = Interpolation(cfg, src_grid, tgt_grid, cache_2);
            Log::info() << interpolation << std::endl;
            interpolation.execute(src_field, tgt_field);
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
    auto improvement = [](double& e, double& r) { return 100. * (r - e) / r; };
    double err;
    // check polygon intersections
    err = remap_stat_1.errors[Statistics::Errors::GEO_DIFF];
    Log::info() << "Polygon area computation improvement: " << improvement(err, tol[0]) << " %" << std::endl;
    EXPECT(err < tol[0]);
    err = remap_stat_1.errors[Statistics::Errors::GEO_L1];
    Log::info() << "Polygon intersection improvement    : " << improvement(err, tol[1]) << " %" << std::endl;
    EXPECT(err < tol[1]);

    // check remap accuracy
    err = remap_stat_1.errors[Statistics::Errors::REMAP_L2];
    Log::info() << "1st order accuracy improvement      : " << improvement(err, tol[2]) << " %" << std::endl;
    EXPECT(err < tol[2]);
    err = remap_stat_2.errors[Statistics::Errors::REMAP_L2];
    Log::info() << "2nd order accuracy improvement      : " << improvement(err, tol[3]) << " %" << std::endl;
    EXPECT(err < tol[3]);

    // check mass conservation
    err = remap_stat_1.errors[Statistics::Errors::REMAP_CONS];
    Log::info() << "1st order conservation improvement  : " << improvement(err, tol[4]) << " %" << std::endl;
    EXPECT(err < tol[4]);
    err = remap_stat_2.errors[Statistics::Errors::REMAP_CONS];
    Log::info() << "2nd order conservation improvement  : " << improvement(err, tol[5]) << " %" << std::endl
                << std::endl;
    EXPECT(err < tol[5]);
    Log::info().unindent();
}

CASE("test_interpolation_conservative") {
#if 1
    SECTION("analytic constfunc") {
        auto func = [](const PointLonLat& p) { return 1.; };
        Statistics remap_stat_1;
        Statistics remap_stat_2;
        do_remapping_test(Grid("H47"), Grid("H48"), func, remap_stat_1, remap_stat_2);
        check(remap_stat_1, remap_stat_2, {1.e-13, 5.e-8, 2.9e-6, 2.9e-6, 5.5e-5, 5.5e-5});
    }

    SECTION("analytic Y_2^2 as in Jones(1998)") {
        auto func = [](const PointLonLat& p) {
            double cos = std::cos(0.025 * p[0]);
            return 2. + cos * cos * std::cos(2 * 0.025 * p[1]);
        };
        Statistics remap_stat_1;
        Statistics remap_stat_2;
        do_remapping_test(Grid("H47"), Grid("H48"), func, remap_stat_1, remap_stat_2);
        check(remap_stat_1, remap_stat_2, {1.e-13, 5.e-8, 4.8e-4, 1.1e-4, 8.9e-5, 1.1e-4});
    }
#endif
}

}  // namespace test
}  // namespace atlas


int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
