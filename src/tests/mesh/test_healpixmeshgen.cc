/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/grid.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/meshgenerator.h"
#include "atlas/meshgenerator/detail/HealpixMeshGenerator.h"
#include "atlas/output/Gmsh.h"
#include "atlas/util/PolygonXY.h"

#include "tests/AtlasTestEnvironment.h"

using namespace atlas::output;
using namespace atlas::meshgenerator;
using namespace atlas::grid;

namespace atlas {
namespace test {

//-----------------------------------------------------------------------------

std::vector<PointLonLat> to_points(std::vector<double>& lon, std::vector<double>& lat) {
    ATLAS_ASSERT(lon.size() == lat.size());
    std::vector<PointLonLat> points(lon.size());
    for (size_t j = 0; j < points.size(); ++j) {
        points[j] = {lon[j], lat[j]};
    }
    return points;
}


CASE("test_check_healpix_points") {
    SECTION("H1") {
        Grid grid("H1");
        EXPECT_EQ(grid.size(), 12);
        std::vector<double> lat = {41.810315, 41.810315, 41.810315,  41.810315,  0.,         0.,
                                   0.,        0.,        -41.810315, -41.810315, -41.810315, -41.810315};
        std::vector<double> lon = {45., 135., 225., 315., 0., 90., 180., 270., 45., 135., 225., 315.};

        auto points = to_points(lon, lat);

        int i = 0;
        for (auto& p : grid.lonlat()) {
            EXPECT_APPROX_EQ(p, points[i], 1.e-6);
            i++;
        }
    }
    SECTION("H2") {
        Grid grid("H2");
        EXPECT_EQ(grid.size(), 48);
        std::vector<double> lat = {66.443536,  66.443536,  66.443536,  66.443536,  41.810315,  41.810315,  41.810315,
                                   41.810315,  41.810315,  41.810315,  41.810315,  41.810315,  19.471221,  19.471221,
                                   19.471221,  19.471221,  19.471221,  19.471221,  19.471221,  19.471221,  0.,
                                   0.,         0.,         0.,         0.,         0.,         0.,         0.,
                                   -19.471221, -19.471221, -19.471221, -19.471221, -19.471221, -19.471221, -19.471221,
                                   -19.471221, -41.810315, -41.810315, -41.810315, -41.810315, -41.810315, -41.810315,
                                   -41.810315, -41.810315, -66.443536, -66.443536, -66.443536, -66.443536};
        std::vector<double> lon = {45.,   135.,  225.,  315.,  22.5,  67.5,  112.5, 157.5, 202.5, 247.5, 292.5, 337.5,
                                   0.,    45.,   90.,   135.,  180.,  225.,  270.,  315.,  22.5,  67.5,  112.5, 157.5,
                                   202.5, 247.5, 292.5, 337.5, 0.,    45.,   90.,   135.,  180.,  225.,  270.,  315.,
                                   22.5,  67.5,  112.5, 157.5, 202.5, 247.5, 292.5, 337.5, 45.,   135.,  225.,  315.};
        auto points             = to_points(lon, lat);

        int i = 0;
        for (auto& p : grid.lonlat()) {
            EXPECT_APPROX_EQ(p, points[i], 1.e-6);
            i++;
        }
    }

    SECTION("H3") {
        Grid grid("H3");
        std::vector<double> lat = {
            74.357529,  74.357529,  74.357529,  74.357529,  58.413662,  58.413662,  58.413662,  58.413662,  58.413662,
            58.413662,  58.413662,  58.413662,  41.810315,  41.810315,  41.810315,  41.810315,  41.810315,  41.810315,
            41.810315,  41.810315,  41.810315,  41.810315,  41.810315,  41.810315,  26.3878,    26.3878,    26.3878,
            26.3878,    26.3878,    26.3878,    26.3878,    26.3878,    26.3878,    26.3878,    26.3878,    26.3878,
            12.839588,  12.839588,  12.839588,  12.839588,  12.839588,  12.839588,  12.839588,  12.839588,  12.839588,
            12.839588,  12.839588,  12.839588,  0.,         0.,         0.,         0.,         0.,         0.,
            0.,         0.,         0.,         0.,         0.,         0.,         -12.839588, -12.839588, -12.839588,
            -12.839588, -12.839588, -12.839588, -12.839588, -12.839588, -12.839588, -12.839588, -12.839588, -12.839588,
            -26.3878,   -26.3878,   -26.3878,   -26.3878,   -26.3878,   -26.3878,   -26.3878,   -26.3878,   -26.3878,
            -26.3878,   -26.3878,   -26.3878,   -41.810315, -41.810315, -41.810315, -41.810315, -41.810315, -41.810315,
            -41.810315, -41.810315, -41.810315, -41.810315, -41.810315, -41.810315, -58.413662, -58.413662, -58.413662,
            -58.413662, -58.413662, -58.413662, -58.413662, -58.413662, -74.357529, -74.357529, -74.357529, -74.357529};

        std::vector<double> lon = {
            45.,  135., 225.,  315.,  22.5,  67.5,  112.5, 157.5, 202.5, 247.5, 292.5, 337.5, 15.,  45.,  75.,  105.,
            135., 165., 195.,  225.,  255.,  285.,  315.,  345.,  0.,    30.,   60.,   90.,   120., 150., 180., 210.,
            240., 270., 300.,  330.,  15.,   45.,   75.,   105.,  135.,  165.,  195.,  225.,  255., 285., 315., 345.,
            0.,   30.,  60.,   90.,   120.,  150.,  180.,  210.,  240.,  270.,  300.,  330.,  15.,  45.,  75.,  105.,
            135., 165., 195.,  225.,  255.,  285.,  315.,  345.,  0.,    30.,   60.,   90.,   120., 150., 180., 210.,
            240., 270., 300.,  330.,  15.,   45.,   75.,   105.,  135.,  165.,  195.,  225.,  255., 285., 315., 345.,
            22.5, 67.5, 112.5, 157.5, 202.5, 247.5, 292.5, 337.5, 45.,   135.,  225.,  315.};

        auto points = to_points(lon, lat);

        int i = 0;
        for (auto& p : grid.lonlat()) {
            EXPECT_APPROX_EQ(p, points[i], 1.e-6);
            i++;
        }
    }
}

CASE("test_create_healpix_mesh") {
    util::Config opts;
    HealpixMeshGenerator generate(opts);

    Mesh m = generate(Grid("H16"));

    Grid grid = m.grid();
    EXPECT(grid);
    Log::info() << grid.spec() << "\n";
    Log::info() << grid.size() << "\n";

    Gmsh("out_3d.msh", util::Config("coordinates", "xyz")).write(m);
    Gmsh("out_ll.msh", util::Config("coordinates", "lonlat")).write(m);
}

//-----------------------------------------------------------------------------

CASE("construction by config") {
    EXPECT(Grid(util::Config("type", "healpix")("N", 3)) == Grid("H3"));
}

CASE("construction by integer") {
    EXPECT(HealpixGrid(3) == Grid("H3"));
}

CASE("construction by integer and ordering") {
    EXPECT(HealpixGrid(3, "ring") == Grid("H3"));
}

CASE("failing construction by integer and unsupported ordering") {
    EXPECT_THROWS(HealpixGrid grid(3, "nested"));
}

CASE("failing construction by config with unsupported ordering") {
    EXPECT_THROWS(Grid grid(util::Config("name", "H3")("ordering", "nested")));
}

CASE("failing construction by config with unsupported ordering") {
    EXPECT_THROWS(Grid grid(util::Config("type", "healpix")("N", 3)("ordering", "nested")));
}

//-----------------------------------------------------------------------------

CASE("matching mesh partitioner") {
    auto grid      = Grid{"H8"};
    auto mesh      = HealpixMeshGenerator{}.generate(grid);
    auto match     = MatchingPartitioner{mesh};
    auto& polygons = mesh.polygons();

    static bool do_once = [&]() {
        for (idx_t i = 0; i < polygons.size(); ++i) {
            auto poly = util::PolygonXY{polygons[i]};
            Log::info() << "polygon[" << i << "]:\n";
            for (idx_t j = 0; j < poly.size(); ++j) {
                Log::info() << "  " << std::setw(5) << std::left << j << std::setprecision(16) << poly[j] << std::endl;
            }
        }
        return true;
    }();


    SECTION("H8 -> O64") { EXPECT_NO_THROW(match.partition(Grid{"O64"})); }
    SECTION("H8 -> L32x17") { EXPECT_NO_THROW(match.partition(Grid{"L32x17"})); }
    SECTION("H8 -> S32x17") { EXPECT_NO_THROW(match.partition(Grid{"S32x17"})); }
    SECTION("H8 -> F32") { EXPECT_NO_THROW(match.partition(Grid{"F32"})); }
    SECTION("H8 -> L64x33 (west=-180)") { EXPECT_NO_THROW(match.partition(Grid{"L64x33", GlobalDomain(-180.)})); }
}

}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
