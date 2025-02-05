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

#include "eckit/types/FloatCompare.h"

#include "atlas/array.h"
#include "atlas/functionspace.h"
#include "atlas/functionspace/PointCloud.h"
#include "atlas/grid.h"
#include "atlas/interpolation.h"
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
//
CASE("test_interpolation_O64_to_empty_PointCloud") {
    Grid grid("O64");
    Mesh mesh(grid);
    NodeColumns fs(mesh);

    atlas::Field points("lonlat", atlas::array::make_datatype<double>(), atlas::array::make_shape(0, 2));
    // Some points at the equator
    PointCloud pointcloud(points);

    Interpolation interpolation(option::type("unstructured-bilinear-lonlat"), fs, pointcloud);

    Field field_source = fs.createField<double>(option::name("source"));
    Field field_target("target", array::make_datatype<double>(), array::make_shape(pointcloud.size()));

    auto source = array::make_view<double, 1>(field_source);
    for (idx_t j = 0; j < fs.nodes().size(); ++j) {
        source(j) = 0;
    }

    interpolation.execute(field_source, field_target);
}

CASE("test_interpolation_O64_to_points_bilinear_remapping") {
    Grid grid("O64");
    Mesh mesh(grid);
    NodeColumns fs(mesh);

    // Some points at the equator
    PointCloud pointcloud(
        {{00., 0.}, {10., 0.}, {20., 0.}, {30., 0.}, {40., 0.}, {50., 0.}, {60., 0.}, {70., 0.}, {80., 0.}, {90., 0.}});

    auto func = [](double x) -> double { return std::sin(x * M_PI / 180.); };

    Interpolation interpolation(
        option::type("unstructured-bilinear-lonlat") | util::Config("max_fraction_elems_to_try", 0.4), fs, pointcloud);

    SECTION("test maximum nearest neighbour settings") {
        std::stringstream test_stream;
        interpolation.print(test_stream);
        std::string test_string = test_stream.str();
        EXPECT((test_string.find("max_fraction_elems_to_try: 0.4") != std::string::npos));
    }

    SECTION("test interpolation outputs") {
        Field field_source = fs.createField<double>(option::name("source"));
        Field field_target("target", array::make_datatype<double>(), array::make_shape(pointcloud.size()));

        auto lonlat = array::make_view<double, 2>(fs.nodes().lonlat());
        auto source = array::make_view<double, 1>(field_source);
        for (idx_t j = 0; j < fs.nodes().size(); ++j) {
            source(j) = func(lonlat(j, LON));
        }

        interpolation.execute(field_source, field_target);

        auto target = array::make_view<double, 1>(field_target);

        auto check = std::vector<double>{func(00.), func(10.), func(20.), func(30.), func(40.),
                                         func(50.), func(60.), func(70.), func(80.), func(90.)};

        for (idx_t j = 0; j < pointcloud.size(); ++j) {
            static double interpolation_tolerance = 1.e-4;
            Log::info() << target(j) << "  " << check[j] << std::endl;
            EXPECT_APPROX_EQ(target(j), check[j], interpolation_tolerance);
        }
    }
}

CASE("test_interpolation_N64_to_O32_bilinear_remapping") {
    Grid grid_src("N64");
    Mesh mesh_src(grid_src);
    NodeColumns fs_src(mesh_src);

    Grid grid_tgt("O32");
    Mesh mesh_tgt(grid_tgt);
    NodeColumns fs_tgt(mesh_tgt);

    Interpolation interpolation(
        option::type("unstructured-bilinear-lonlat") | util::Config("max_fraction_elems_to_try", 0.4), fs_src, fs_tgt);

    Field field_source = fs_src.createField<double>(option::name("source"));
    Field field_target = fs_tgt.createField<double>(option::name("target"));

    const double deg2rad = M_PI / 180.;
    auto func = [](double lon, double lat, double t) { return std::cos(lat) * std::sin(lon); };

    array::ArrayView<double, 2> lonlat = array::make_view<double, 2>(fs_src.nodes().lonlat());
    array::ArrayView<double, 1> source = array::make_view<double, 1>(field_source);

    for (idx_t j = 0; j < lonlat.shape(0); ++j) {
        const double lon = deg2rad * lonlat(j, 0);
        const double lat = deg2rad * lonlat(j, 1);
        source(j)        = func(lon, lat, 1.);
    }

    interpolation.execute(field_source, field_target);

    array::ArrayView<double, 1> target     = array::make_view<double, 1>(field_target);
    array::ArrayView<double, 2> lonlat_tgt = array::make_view<double, 2>(fs_tgt.nodes().lonlat());
    std::vector<double> diffs;
    ATLAS_ASSERT(target.shape(0) == lonlat_tgt.shape(0));
    for (idx_t j = 0; j < lonlat_tgt.shape(0); ++j) {
        const double lon   = deg2rad * lonlat_tgt(j, 0);
        const double lat   = deg2rad * lonlat_tgt(j, 1);
        const double check = func(lon, lat, 1.);
        diffs.push_back((target(j) - check) * (target(j) - check));
        if (std::abs(target(j) - check) > 0.1) {
            Log::info() << "lon, lat:" << lon << ", " << lat << " target, check: " << target(j) << ", " << check
                        << std::endl;
        }
    }
    static double interpolation_tolerance = 1.e-3;
    double std_dev = std::accumulate(diffs.begin(), diffs.end(), decltype(diffs)::value_type(0)) / diffs.size();
    double max_dev = std::sqrt(*std::max_element(diffs.begin(), diffs.end()));
    Log::info() << " standard deviation " << std_dev << " max " << max_dev << std::endl;
    EXPECT_APPROX_EQ(std_dev, 0.0, 1.e-6);
    EXPECT_APPROX_EQ(max_dev, 0.0, interpolation_tolerance);
}
//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
