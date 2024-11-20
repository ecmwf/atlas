/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <cmath>

#include "eckit/log/Bytes.h"

#include "atlas/array.h"
#include "atlas/functionspace.h"
#include "atlas/functionspace/PointCloud.h"
#include "atlas/grid.h"
#include "atlas/interpolation.h"
#include "atlas/linalg/sparse.h"
#include "atlas/mesh.h"
#include "atlas/meshgenerator.h"
#include "atlas/util/CoordinateEnums.h"

#include "tests/AtlasTestEnvironment.h"

using namespace eckit;
using namespace atlas::functionspace;
using namespace atlas::util;

namespace atlas {
namespace test {

//-----------------------------------------------------------------------------

void set_field(Field& field, Grid& grid, std::function<double(double, double)> func) {
    auto view = array::make_view<double, 1>(field);
    idx_t j{0};
    for (auto& p : grid.lonlat()) {
        view(j) = func(p.lon(), p.lat());
        ++j;
    }
}

void set_field(Field& field, double v) {
    auto view = array::make_view<double, 1>(field);
    for (idx_t j = 0; j < view.shape(0); ++j) {
        view(j) = v;
    }
}


void check_field(Field& field, Grid& grid, std::function<double(double, double)> func, double tolerance) {
    auto view = array::make_view<double, 1>(field);
    idx_t j{0};
    for (auto& p : grid.lonlat()) {
        EXPECT_APPROX_EQ(view(j), func(p.lon(), p.lat()), tolerance);
        ++j;
    }
}

interpolation::MatrixCache get_or_create_cache(Grid& source, Grid& target) {
    ATLAS_TRACE();
    static std::map<std::string, interpolation::Cache> map;
    std::string key = source.uid() + target.uid();
    if (map.find(key) == map.end()) {
        // Create interpolator, extract matrix, and store it in map
        map[key] = interpolation::MatrixCache(Interpolation(option::type("finite-element"), source, target));
    }
    return map[key];
}

auto func = [](double x, double) -> double { return std::sin(x * M_PI / 180.); };

//-----------------------------------------------------------------------------

CASE("extract cache and use") {
    Grid grid_source("F32");
    Grid grid_target("F16");

    Field field_source("source", array::make_datatype<double>(), array::make_shape(grid_source.size()));
    Field field_target("target", array::make_datatype<double>(), array::make_shape(grid_target.size()));

    set_field(field_source, grid_source, func);

    interpolation::Cache cache = get_or_create_cache(grid_source, grid_target);

    Log::info() << "Cache contains " << eckit::Bytes(cache.footprint()) << std::endl;

    ATLAS_TRACE_SCOPE("Interpolate with cache") {
        Interpolation interpolation_using_cache(option::type("finite-element"), grid_source, grid_target, cache);
        interpolation_using_cache.execute(field_source, field_target);
    }

    check_field(field_target, grid_target, func, 1.e-4);
}

//-----------------------------------------------------------------------------

CASE("extract cache, copy it, and move it for use") {
    Grid grid_source("F32");
    Grid grid_target("F16");

    Field field_source("source", array::make_datatype<double>(), array::make_shape(grid_source.size()));
    Field field_target("target", array::make_datatype<double>(), array::make_shape(grid_target.size()));

    set_field(field_source, grid_source, func);

    atlas::linalg::SparseMatrix matrix = get_or_create_cache(grid_source, grid_target).matrix();

    EXPECT(not matrix.empty());

    auto cache = interpolation::MatrixCache(std::move(matrix));

    EXPECT(matrix.empty());

    ATLAS_TRACE_SCOPE("Interpolate with cache") {
        Interpolation interpolation_using_cache(option::type("finite-element"), grid_source, grid_target, cache);
        interpolation_using_cache.execute(field_source, field_target);
    }

    check_field(field_target, grid_target, func, 1.e-4);
}

//-----------------------------------------------------------------------------

CASE("extract cache, copy it, and pass non-owning pointer") {
    Grid grid_source("F32");
    Grid grid_target("F16");

    Field field_source("source", array::make_datatype<double>(), array::make_shape(grid_source.size()));
    Field field_target("target", array::make_datatype<double>(), array::make_shape(grid_target.size()));

    set_field(field_source, grid_source, func);

    atlas::linalg::SparseMatrix matrix = get_or_create_cache(grid_source, grid_target).matrix();

    EXPECT(not matrix.empty());

    ATLAS_TRACE_SCOPE("Interpolate with [eckit_linalg] sparse_matrix_multiply of ArrayView") {
        atlas::linalg::sparse::current_backend("eckit_linalg");
        atlas::linalg::sparse::default_backend("eckit_linalg").set("backend", "generic");
        auto src = array::make_view<double, 1>(field_source);
        auto tgt = array::make_view<double, 1>(field_target);
        atlas::linalg::sparse_matrix_multiply(matrix, src, tgt);
    }
    check_field(field_target, grid_target, func, 1.e-4);
    set_field(field_target, 0.);
    ATLAS_TRACE_SCOPE("Interpolate with [openmp] sparse_matrix_multiply of eckit::linalg::Vector") {
        atlas::linalg::sparse::current_backend("openmp");
        eckit::linalg::Vector src{array::make_view<double, 1>(field_source).data(), field_source.size()};
        eckit::linalg::Vector tgt{array::make_view<double, 1>(field_target).data(), field_target.size()};
        atlas::linalg::sparse_matrix_multiply(matrix, src, tgt);
    }
    check_field(field_target, grid_target, func, 1.e-4);
    set_field(field_target, 0.);

    auto cache = interpolation::MatrixCache(&matrix);

    EXPECT(not matrix.empty());

    ATLAS_TRACE_SCOPE("Interpolate with cache") {
        set_field(field_target, 0.);
        Interpolation interpolation_using_cache(option::type("finite-element"), grid_source, grid_target, cache);
        interpolation_using_cache.execute(field_source, field_target);
    }
    check_field(field_target, grid_target, func, 1.e-4);
    set_field(field_target, 0.);
}

//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
