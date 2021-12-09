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

CASE("test_interpolation_finite_element") {
    Grid grid("O64");
    Mesh mesh(grid);
    NodeColumns fs(mesh);

    // Some points at the equator
    PointCloud pointcloud(
        {{00., 0.}, {10., 0.}, {20., 0.}, {30., 0.}, {40., 0.}, {50., 0.}, {60., 0.}, {70., 0.}, {80., 0.}, {90., 0.}});

    auto func = [](double x) -> double { return std::sin(x * M_PI / 180.); };

    Interpolation interpolation(option::type("finite-element") | util::Config("max_fraction_elems_to_try", 0.4), fs,
                                pointcloud);

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
            EXPECT(eckit::types::is_approximately_equal(target(j), check[j], interpolation_tolerance));
        }
    }
}

//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
