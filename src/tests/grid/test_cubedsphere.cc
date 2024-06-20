/*
 * (C) Crown Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "atlas/array/MakeView.h"
#include "atlas/field/FieldSet.h"
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/grid.h"
#include "atlas/grid/Partitioner.h"
#include "atlas/grid/Tiles.h"
#include "atlas/grid/detail/partitioner/CubedSpherePartitioner.h"
#include "atlas/grid/detail/tiles/Tiles.h"
#include "atlas/mesh.h"
#include "atlas/meshgenerator.h"
#include "atlas/option.h"
#include "atlas/output/Gmsh.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/util/CoordinateEnums.h"

#include "tests/AtlasTestEnvironment.h"

namespace atlas {
namespace test {
namespace {

using grid::detail::partitioner::CubedSpherePartitioner;

void partition(const CubedSpherePartitioner& partitioner, const Grid& grid, CubedSpherePartitioner::CubedSphere& cb,
               std::vector<int>& part) {
    std::vector<CubedSpherePartitioner::CellInt> nodes(static_cast<std::size_t>(grid.size()));
    std::size_t n(0);

    for (std::size_t it = 0; it < 6; ++it) {
        for (idx_t iy = 0; iy < cb.ny[it]; ++iy) {
            for (idx_t ix = 0; ix < cb.nx[it]; ++ix) {
                nodes[n].t = static_cast<int>(it);
                nodes[n].x = static_cast<int>(ix);
                nodes[n].y = static_cast<int>(iy);
                nodes[n].n = static_cast<int>(n);
                ++n;
            }
        }
    }
    partitioner.partition(cb, grid.size(), nodes.data(), part.data());
}
}  // namespace

CASE("cubedsphere_tile_test") {
    auto tileConfig1 = atlas::util::Config("type", "cubedsphere_lfric");
    auto lfricTiles  = atlas::grid::CubedSphereTiles(tileConfig1);
    EXPECT(lfricTiles.type() == "cubedsphere_lfric");

    auto tileConfig2 = atlas::util::Config("type", "cubedsphere_fv3");
    auto fv3Tiles    = atlas::grid::CubedSphereTiles(tileConfig2);
    EXPECT(fv3Tiles.type() == "cubedsphere_fv3");

    auto lfricTiles2 = atlas::grid::CubedSphereTiles("cubedsphere_lfric");
    EXPECT(lfricTiles2.type() == "cubedsphere_lfric");

    auto fv3Tiles2 = atlas::grid::CubedSphereTiles("cubedsphere_fv3");
    EXPECT(fv3Tiles.type() == "cubedsphere_fv3");
}


//-----------------------------------------------------------------------------

CASE("test_iterator") {
    std::vector<int> resolutions{1, 2, 4, 8};
    std::vector<std::string> grid_prefixes{"CS-EA-L-", "CS-ED-L-", "CS-EA-C-", "CS-ED-C-", "CS-LFR-L-", "CS-LFR-C-"};


    for (auto resolution : resolutions) {
        for (auto& grid_prefix : grid_prefixes) {
            std::string grid_name = grid_prefix + std::to_string(resolution);
            SECTION(grid_name) {
                if (grid_name == "CS-LFR-L-1") {
                    Log::error() << eckit::Colour::red << "TODO: Fix me!!!. Skipping..." << eckit::Colour::reset
                                 << std::endl;
                    continue;
                }
                Grid g(grid_name);

                // Test xy
                {
                    std::vector<PointXY> coordinates_1;
                    std::vector<PointXY> coordinates_2;
                    std::vector<PointXY> coordinates_3;
                    {
                        for (auto crd : g.xy()) {
                            coordinates_1.push_back(crd);
                        }
                    }
                    {
                        auto iterator = g.xy().begin();
                        for (int n = 0; n < g.size(); ++n) {
                            coordinates_2.push_back(*iterator);
                            iterator += 1;
                        }
                    }
                    {
                        auto iterator = g.xy().begin();
                        PointXY crd;
                        while (iterator.next(crd)) {
                            coordinates_3.push_back(crd);
                        }
                    }
                    EXPECT_EQ(coordinates_1.size(), g.size());
                    EXPECT_EQ(coordinates_2.size(), g.size());
                    EXPECT_EQ(coordinates_3.size(), g.size());
                    EXPECT_EQ(coordinates_2, coordinates_1);
                    EXPECT_EQ(coordinates_3, coordinates_1);
                }

                // Test lonlat
                {
                    std::vector<PointLonLat> coordinates_1;
                    std::vector<PointLonLat> coordinates_2;
                    std::vector<PointLonLat> coordinates_3;
                    {
                        for (auto crd : g.lonlat()) {
                            coordinates_1.push_back(crd);
                        }
                    }
                    {
                        auto iterator = g.lonlat().begin();
                        for (int n = 0; n < g.size(); ++n) {
                            coordinates_2.push_back(*iterator);
                            iterator += 1;
                        }
                    }
                    {
                        auto iterator = g.lonlat().begin();
                        PointLonLat crd;
                        while (iterator.next(crd)) {
                            coordinates_3.push_back(crd);
                        }
                    }
                    EXPECT_EQ(coordinates_1.size(), g.size());
                    EXPECT_EQ(coordinates_2.size(), g.size());
                    EXPECT_EQ(coordinates_3.size(), g.size());
                    EXPECT_EQ(coordinates_2, coordinates_1);
                    EXPECT_EQ(coordinates_3, coordinates_1);
                }
            }
        }
    }
}

CASE("cubedsphere_FV3_mesh_test") {
    // THIS IS TEMPORARY!
    // I expect this will be replaced by some more aggressive tests.

    // Set grid.
    const auto grid = atlas::Grid("CS-EA-L-2");

    atlas::Log::info() << grid.name() << std::endl;
    atlas::Log::info() << grid.size() << std::endl;


    // Set mesh.
    auto meshGen = atlas::MeshGenerator("nodal-cubedsphere");
    auto mesh    = meshGen.generate(grid);

    // Set functionspace
    auto functionSpace = atlas::functionspace::NodeColumns(mesh);

    // Set field
    auto field = functionSpace.ghost();


    // Set gmsh config.
    auto gmshConfigXy     = atlas::util::Config("coordinates", "xy") | atlas::util::Config("ghost", false);
    auto gmshConfigXyz    = atlas::util::Config("coordinates", "xyz") | atlas::util::Config("ghost", false);
    auto gmshConfigLonLat = atlas::util::Config("coordinates", "lonlat") | atlas::util::Config("ghost", false);

    // Set source gmsh object.
    const auto gmshXy     = atlas::output::Gmsh("FV3_xy_mesh.msh", gmshConfigXy);
    const auto gmshXyz    = atlas::output::Gmsh("FV3_xyz_mesh.msh", gmshConfigXyz);
    const auto gmshLonLat = atlas::output::Gmsh("FV3_lonlat_mesh.msh", gmshConfigLonLat);

    // Write gmsh.
    gmshXy.write(mesh);
    gmshXy.write(field);
    gmshXyz.write(mesh);
    gmshXyz.write(field);
    gmshLonLat.write(mesh);
    gmshLonLat.write(field);
}

CASE("cubedsphere_generic_mesh_test") {
    // THIS IS TEMPORARY!
    // I expect this will be replaced by some more aggressive tests.

    // Set grid.
    const auto grid = atlas::Grid("CS-LFR-C-32");

    atlas::Log::info() << grid.name() << std::endl;
    atlas::Log::info() << grid.size() << std::endl;


    // Set mesh.
    auto meshGen = atlas::MeshGenerator("cubedsphere");
    auto mesh    = meshGen.generate(grid);


    // Visually inspect the fields.
    // remote indices should be equal at tile boundaries.


    // Set gmsh config.
    auto gmshConfigXy     = atlas::util::Config("coordinates", "xy");
    auto gmshConfigXyz    = atlas::util::Config("coordinates", "xyz");
    auto gmshConfigLonLat = atlas::util::Config("coordinates", "lonlat");

    gmshConfigXy.set("ghost", true);
    gmshConfigXy.set("info", true);

    gmshConfigXyz.set("ghost", true);
    gmshConfigXyz.set("info", true);

    gmshConfigLonLat.set("ghost", true);
    gmshConfigLonLat.set("info", true);

    // Set source gmsh object.
    const auto gmshXy     = atlas::output::Gmsh("cs_xy_mesh.msh", gmshConfigXy);
    const auto gmshXyz    = atlas::output::Gmsh("cs_xyz_mesh.msh", gmshConfigXyz);
    const auto gmshLonLat = atlas::output::Gmsh("cs_lonlat_mesh.msh", gmshConfigLonLat);

    // Write gmsh.
    gmshXy.write(mesh);
    gmshXyz.write(mesh);
    gmshLonLat.write(mesh);
}

CASE("cubedsphere_tileCubePeriodicity_test") {
    auto tileConfig1 = atlas::util::Config("type", "cubedsphere_lfric");
    auto lfricTiles  = atlas::grid::CubedSphereTiles(tileConfig1);

    // create a nodal cubed-sphere grid and check that no point are changed by
    // iterating through points.

    int resolution(2);
    std::vector<std::string> grid_names{
        "CS-LFR-L-" + std::to_string(resolution),
    };
    Grid grid{grid_names[0]};

    int jn{0};
    for (auto crd : grid.xy()) {
        atlas::PointXY initialXY{crd[XX], crd[YY]};
        double xy[2]           = {initialXY.x(), initialXY.y()};
        idx_t t                = lfricTiles.indexFromXY(xy);
        atlas::PointXY finalXY = lfricTiles.tileCubePeriodicity(initialXY, t);
        EXPECT_APPROX_EQ(initialXY, finalXY);
        ++jn;
    }

    std::vector<atlas::PointXY> startingXYTile0{{0., 315.},  {90., 315.},  {0., 225.},   {90., 225.},  {0., 135.},
                                                {90., 135.}, {45., 0.},    {45., 90.},   {0., -45.},   {90., -45.},
                                                {0, -135.},  {90, -135.},  {-90., 45.},  {180., 45.},  {270., 45.},
                                                {360., 45.}, {-90., -45.}, {180., -45.}, {270., -45.}, {360., -45.}};

    std::vector<atlas::PointXY> expectedXYTile0{{0., -45.},   {90., -45.},  {270., -45.}, {180., -45.}, {0., 135.},
                                                {90., 135.},  {45., 0.},    {45., 90.},   {0., -45.},   {90., -45.},
                                                {270., -45.}, {180., -45.}, {0., 135.},   {90., 135.},  {0., 135.},
                                                {0., 45.},    {270., -45.}, {180., -45.}, {270., -45.}, {0., -45.}};

    std::vector<atlas::PointXY> startingXYTile1{{90., 315.},  {180., 315.}, {90., 225.},  {180., 225.}, {90., 135.},
                                                {180., 135.}, {135., 0.},   {135., 90.},  {90., -45.},  {180., -45.},
                                                {90., -135.}, {180, -135.}, {0., 45.},    {270., 45.},  {360., 45.},
                                                {450., 45.},  {0., -45.},   {270., -45.}, {360., -45.}, {450., -45.}};

    std::vector<atlas::PointXY> expectedXYTile1{{90., -45.}, {180., -45.}, {0., -45.},   {270., -45.}, {0., 45.},
                                                {0., 135.},  {135., 0.},   {45., 90.},   {90., -45.},  {180., -45.},
                                                {0., -45.},  {270., -45.}, {0., 45.},    {0., 135.},   {0., 45.},
                                                {90., 45.},  {0., -45.},   {270., -45.}, {0., -45.},   {90., -45.}};


    std::vector<atlas::PointXY> startingXYTile2{{180., 315.},  {270., 315.}, {180., 225.}, {270., 225.}, {180., 135.},
                                                {270., 135.},  {225., 0.},   {225., 90.},  {180., -45.}, {270., -45.},
                                                {180., -135.}, {270, -135.}, {90., 45.},   {0., 45.},    {-90., 45.},
                                                {180., 45.},   {90., -45.},  {0., -45.},   {90., -45.},  {180., -45.}};

    std::vector<atlas::PointXY> expectedXYTile2{{180., -45.},  {270., -45.}, {90., -45.}, {0., -45.},   {90., 45.},
                                                {0., 45.},     {225., 0.},   {45., 90.},  {180., -45.}, {270., -45.},
                                                {90., -45.},   {0., -45.},   {90., 45.},  {0., 45.},    {0., 135.},
                                                {90.0, 135.0}, {90., -45.},  {0., -45.},  {90., -45.},  {180., -45.}};

    std::vector<atlas::PointXY> startingXYTile3{{270., 315.},  {360., 315.}, {270., 225.}, {360., 225.}, {270., 135.},
                                                {360., 135.},  {315., 90.},  {315., -90.}, {270., -45.}, {360., -45.},
                                                {270., -135.}, {360, -135.}, {180., 45.},  {90., 45.},   {360., 45.},
                                                {270., 45.},   {180., -45.}, {90., -45.},  {360., -45.}, {270., -45.}};

    std::vector<atlas::PointXY> expectedXYTile3{{270., -45.}, {0., -45.},   {180., -45.}, {90., -45.},  {90., 135.},
                                                {90., 45.},   {45., 90.},   {45., -90.},  {270., -45.}, {0., -45.},
                                                {180., -45.}, {90., -45.},  {90., 135.},  {90., 45.},   {0., 45.},
                                                {0.0, 135.0}, {180., -45.}, {90., -45.},  {0., -45.},   {270., -45.}};


    std::vector<atlas::PointXY> startingXYTile4{{0., 405.},   {90., 405.}, {0., 315.},   {90., 315.},  {0., 225.},
                                                {90., 225.},  {0, -135.},  {45., -90.},  {0., 45.},    {90., 45.},
                                                {0, -45.},    {90, -45.},  {-90., 135.}, {180., 135.}, {270., 135.},
                                                {360., 135.}, {-90., 45.}, {180., 45.},  {270., 45.},  {360., 45.}};

    std::vector<atlas::PointXY> expectedXYTile4{{0., 45.},    {90., 45.},   {0., -45.},   {90., -45.},  {270., -45.},
                                                {180., -45.}, {270., -45.}, {45., -90.},  {0., 45.},    {90., 45.},
                                                {0., -45.},   {90., -45.},  {270., -45.}, {180., -45.}, {270., -45.},
                                                {0., 135.},   {0., -45.},   {90., -45.},  {0., -45.},   {0., 45.}};

    std::vector<atlas::PointXY> startingXYTile5{
        {0., 225.},   {90., 225.},  {0., 135.},    {90., 135.},   {0., 45.},     {90., 45.},   {45., -90.},
        {45., 0.},    {0., -135.},  {90., -135.},  {0, -225.},    {90, -225.},   {-90., -45.}, {180., -45.},
        {270., -45.}, {360., -45.}, {-90., -135.}, {180., -135.}, {270., -135.}, {360., -135.}};

    std::vector<atlas::PointXY> expectedXYTile5{{270., -45.}, {180., -45.}, {0., 135.},  {90., 135.},  {0., 45.},
                                                {90., 45.},   {45., -90.},  {45., 0.},   {270., -45.}, {180., -45.},
                                                {0., 135.},   {90., 135.},  {0., 45.},   {90., 45.},   {0., 45.},
                                                {0., -45.},   {0., 135.},   {90., 135.}, {0., 135.},   {270., -45.}};


    // testing tile 0
    for (atlas::idx_t t = 0; t < 6; t++) {
        std::vector<atlas::PointXY> startingXY;
        std::vector<atlas::PointXY> expectedXY;

        if (t == 0) {
            startingXY = startingXYTile0;
            expectedXY = expectedXYTile0;
        }
        if (t == 1) {
            startingXY = startingXYTile1;
            expectedXY = expectedXYTile1;
        }
        if (t == 2) {
            startingXY = startingXYTile2;
            expectedXY = expectedXYTile2;
        }
        if (t == 3) {
            startingXY = startingXYTile3;
            expectedXY = expectedXYTile3;
        }
        if (t == 4) {
            startingXY = startingXYTile4;
            expectedXY = expectedXYTile4;
        }
        if (t == 5) {
            startingXY = startingXYTile5;
            expectedXY = expectedXYTile5;
        }

        std::size_t jn{0};
        for (atlas::PointXY p : startingXY) {
            atlas::PointXY middleXY = lfricTiles.tileCubePeriodicity(p, t);
            atlas::PointXY finalXY  = lfricTiles.tileCubePeriodicity(middleXY, 0);
            EXPECT_APPROX_EQ(middleXY, finalXY);
            EXPECT_APPROX_EQ(middleXY, expectedXY[jn]);
            ++jn;
        }
    }
}

CASE("cubedsphere_partitioner_test") {
    int resolution(4);
    std::vector<std::string> grid_names{
        "CS-LFR-C-" + std::to_string(resolution),
    };
    Grid grid{grid_names[0]};

    using grid::detail::partitioner::CubedSpherePartitioner;

    if (mpi::size() == 1) {
        // factory based constructor
        SECTION("factory based constructor")
        {
            std::vector<int> globalProcStartPE{0, 0, 0, 1, 1, 1};
            std::vector<int> globalProcEndPE{0, 0, 0, 1, 1, 1};
            std::vector<int> nprocx{1, 1, 1, 1, 1, 1};
            std::vector<int> nprocy{1, 1, 1, 1, 1, 1};

            atlas::util::Config conf;
            conf.set("starting rank on tile", globalProcStartPE);
            conf.set("final rank on tile", globalProcEndPE);
            conf.set("nprocx", nprocx);
            conf.set("nprocy", nprocy);

            grid::Partitioner partitioner("cubedsphere", conf);
            grid::Distribution d_cs = partitioner.partition(grid);

            grid::Partitioner partitioner2("cubedsphere", 2);
            grid::Distribution d_cs2 = partitioner2.partition(grid);
        }

        // 2 partitions via configuration and distribution object
        SECTION("2 partitions via configuration and distribution object")
        {
            std::vector<int> globalProcStartPE{0, 0, 0, 1, 1, 1};
            std::vector<int> globalProcEndPE{0, 0, 0, 1, 1, 1};
            std::vector<int> nprocx{1, 1, 1, 1, 1, 1};
            std::vector<int> nprocy{1, 1, 1, 1, 1, 1};

            atlas::util::Config conf;
            conf.set("starting rank on tile", globalProcStartPE);
            conf.set("final rank on tile", globalProcEndPE);
            conf.set("nprocx", nprocx);
            conf.set("nprocy", nprocy);

            grid::Distribution d_cs = grid::Partitioner(new CubedSpherePartitioner(2, conf)).partition(grid);

            for (idx_t t = 0; t < 96; ++t) {
                EXPECT(d_cs.partition(t) == t / 48);
            }
        }

        // 3 partitions  via vector constructor and distribution object
        SECTION("3 partitions")
        {
            std::vector<int> globalProcStartPE{0, 0, 1, 1, 2, 2};
            std::vector<int> globalProcEndPE{0, 0, 1, 1, 2, 2};
            std::vector<int> nprocx{1, 1, 1, 1, 1, 1};
            std::vector<int> nprocy{1, 1, 1, 1, 1, 1};

            grid::Distribution d_cs =
                grid::Partitioner(new CubedSpherePartitioner(3, globalProcStartPE, globalProcEndPE, nprocx, nprocy))
                    .partition(grid);

            for (idx_t t = 0; t < 96; ++t) {
                EXPECT(d_cs.partition(t) == t / 32);
            }
        }

        // 4 partitions
        SECTION("4 partitions")
        {
            CubedSpherePartitioner partitioner(4);
            CubedSpherePartitioner::CubedSphere cb = partitioner.cubedsphere(grid);
            std::vector<int> part(static_cast<size_t>(grid.size()), 0);
            partition(partitioner, grid, cb, part);

            for (std::size_t t = 0; t < 4; ++t) {
                EXPECT(cb.nproc[t] == 1);
                EXPECT(cb.nprocx[t] == 1);
                EXPECT(cb.nprocy[t] == 1);
                EXPECT(cb.nx[t] == 4);
                EXPECT(cb.ny[t] == 4);
                EXPECT(cb.globalProcStartPE[t] == static_cast<atlas::idx_t>(t));
                EXPECT(cb.globalProcEndPE[t] == static_cast<atlas::idx_t>(t));
            }

            for (std::size_t t = 4; t < 6; ++t) {
                EXPECT(cb.nproc[t] == 0);
                EXPECT(cb.nprocx[t] == 1);
                EXPECT(cb.nprocy[t] == 1);
                EXPECT(cb.nx[t] == 4);
                EXPECT(cb.ny[t] == 4);
                EXPECT(cb.globalProcStartPE[t] == static_cast<atlas::idx_t>(3));
                EXPECT(cb.globalProcEndPE[t] == static_cast<atlas::idx_t>(3));
            }

            for (size_t i = 0; i < static_cast<size_t>(grid.size()); ++i) {
                if (i < 64) {
                    EXPECT(part[i] == static_cast<idx_t>(i / 16));
                }
                else {
                    EXPECT(part[i] == 3);
                }
            }
        }

        // 12 partitions
        SECTION("12 partitions")
        {
            CubedSpherePartitioner partitioner(12);
            CubedSpherePartitioner::CubedSphere cb = partitioner.cubedsphere(grid);
            std::vector<int> part(static_cast<size_t>(grid.size()), 0);
            partition(partitioner, grid, cb, part);

            for (std::size_t t = 0; t < 6; ++t) {
                EXPECT(cb.nproc[t] == 2);
                EXPECT(cb.nprocx[t] == 1);
                EXPECT(cb.nprocy[t] == 2);
                EXPECT(cb.nx[t] == 4);
                EXPECT(cb.ny[t] == 4);
                EXPECT(cb.globalProcStartPE[t] == static_cast<int>(2 * t));
                EXPECT(cb.globalProcEndPE[t] == static_cast<int>(2 * t + 1));

                for (size_t i = 0; i < static_cast<size_t>(grid.size()); ++i) {
                    EXPECT(part[i] == static_cast<int>(i / 8));
                }
            }
        }

        // 24 partitions
        SECTION("24 partitions")
        {
            CubedSpherePartitioner partitioner(24);
            CubedSpherePartitioner::CubedSphere cb = partitioner.cubedsphere(grid);
            std::vector<int> part(static_cast<size_t>(grid.size()), 0);
            partition(partitioner, grid, cb, part);

            for (std::size_t t = 0; t < 6; ++t) {
                EXPECT(cb.nproc[t] == 4);
                EXPECT(cb.nprocx[t] == 2);
                EXPECT(cb.nprocy[t] == 2);
                EXPECT(cb.nx[t] == 4);
                EXPECT(cb.ny[t] == 4);
                EXPECT(cb.globalProcStartPE[t] == static_cast<int>(4 * t));
                EXPECT(cb.globalProcEndPE[t] == static_cast<int>(4 * t + 3));

                std::size_t l(0);
                for (idx_t t = 0; t < 6; ++t) {
                    for (idx_t j = 0; j < 4; ++j) {
                        for (idx_t i = 0; i < 4; ++i, ++l) {
                            gidx_t temp  = l / 2;
                            gidx_t temp2 = temp % 2;
                            gidx_t temp8 = l / 8;
                            EXPECT(part[l] == temp2 + 2 * temp8);
                        }
                    }
                }
            }
        }

        // 24 partitions, creating distribution object
        SECTION("24 partitions, creating distribution object")
        {
            grid::Distribution d_cs = grid::Partitioner(new CubedSpherePartitioner(24)).partition(grid);
            gidx_t l(0);
            for (idx_t t = 0; t < 6; ++t) {
                for (idx_t j = 0; j < 4; ++j) {
                    for (idx_t i = 0; i < 4; ++i, ++l) {
                        gidx_t temp  = l / 2;
                        gidx_t temp2 = temp % 2;
                        gidx_t temp8 = l / 8;
                        EXPECT(d_cs.partition(l) == temp2 + 2 * temp8);
                    }
                }
            }
        }
    }
    else {
        // Factory test using mpi to determine number of partitions
        grid::Partitioner partitioner2("cubedsphere", atlas::mpi::size());
        grid::Distribution d_cs2 = partitioner2.partition(grid);
    }
}
}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
