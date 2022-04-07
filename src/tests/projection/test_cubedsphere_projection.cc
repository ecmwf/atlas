/*
 * (C) Crown Copyright 2022 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "atlas/array/MakeView.h"
#include "atlas/field/FieldSet.h"
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/grid.h"
#include "atlas/grid/CubedSphereGrid.h"
#include "atlas/mesh.h"
#include "atlas/meshgenerator.h"
#include "atlas/output/Gmsh.h"
#include "atlas/projection/detail/CubedSphereProjectionBase.h"
#include "atlas/util/Constants.h"
#include "tests/AtlasTestEnvironment.h"
namespace atlas {
namespace test {

// Small number relative to 360.
constexpr double epsilon = std::numeric_limits<double>::epsilon() * 360.;

void testProjection(const std::string& gridType, const std::string& meshType, const std::string& outputID) {
    // Create grid.
    const auto grid = CubedSphereGrid(gridType);

    // Get projection.
    const auto& csProjection = grid.cubedSphereProjection();

    // Output tile centres and Jacobains.
    for (size_t i = 0; i < 6; ++i) {
        Log::info() << "Tile " << i << std::endl
                    << "Centre:" << std::endl
                    << csProjection.getCubedSphereTiles().tileCentre(i) << std::endl
                    << "Jacobian:" << std::endl
                    << csProjection.getCubedSphereTiles().tileJacobian(i) << std::endl
                    << std::endl;
    }

    // Create mesh.
    auto mesh = MeshGenerator(meshType, util::Config("halo", 3)).generate(grid);

    // Output mesh.
    const auto gmshConfig =
        util::Config("coordinates", "xy") | util::Config("ghost", true) | util::Config("info", true);
    auto gmsh = output::Gmsh(outputID + "_before.msh", gmshConfig);
    gmsh.write(mesh);


    // "Hack" mesh xy coordinates.
    auto xyView        = array::make_view<double, 2>(mesh.nodes().xy());
    const auto tijView = array::make_view<idx_t, 2>(mesh.nodes().field("tij"));
    for (idx_t i = 0; i < mesh.nodes().size(); ++i) {
        const auto t = tijView(i, 0);

        const auto xy        = PointXY(xyView(i, XX), xyView(i, YY));
        const auto alphabeta = csProjection.xy2alphabeta(xy, t);

        // Inverse function is degenerate when abs(alpha) == abs(beta) and
        // abs(alpha) > 45.
        const bool degenerate =
            std::abs(alphabeta[0]) > 45. && approx_eq(std::abs(alphabeta[0]), std::abs(alphabeta[1]));

        // Check inverse function.
        if (!degenerate) {
            const auto newXy = csProjection.alphabeta2xy(alphabeta, t);
            EXPECT_APPROX_EQ(xy[XX], newXy[XX], epsilon);
            EXPECT_APPROX_EQ(xy[YY], newXy[YY], epsilon);
        }

        // overwrite mesh xy field.
        xyView(i, XX) = alphabeta[0];
        xyView(i, YY) = alphabeta[1];
    }

    // Output mesh updated mesh.
    gmsh = output::Gmsh(outputID + "_after.msh", gmshConfig);
    gmsh.write(mesh);
}

void testJacobian(const std::string& gridType, const std::string& meshType, const std::string& outputID) {
    // Create grid.
    const auto grid = CubedSphereGrid(gridType);

    // Get projection.
    const auto& csProjection = grid.cubedSphereProjection();

    // Set xy0 values on different tiles.
    const auto xy0 = std::array<PointXY, 6>{PointXY{10., -35},  PointXY{100., -35}, PointXY{190., -35},
                                            PointXY{280., -35}, PointXY{10., 55},   PointXY{10., -125}};

    // Get equivalent lonlat values.
    const auto lonlat0 = std::array<PointLonLat, 6>{csProjection.lonlat(xy0[0]), csProjection.lonlat(xy0[1]),
                                                    csProjection.lonlat(xy0[2]), csProjection.lonlat(xy0[3]),
                                                    csProjection.lonlat(xy0[4]), csProjection.lonlat(xy0[5])};
    // Set lonlat step size.
    const auto hLonlat0 = PointLonLat{5., 10.};

    // Perform first order extrapolation test on each tile.
    for (size_t t = 0; t < 6; ++t) {
        const auto jac = csProjection.jacobian(lonlat0[t]);

        double dxyOld{};

        for (int i = 0; i < 5; ++i) {
            // Sequentially halve lonlat increment.
            const auto hLonlat = hLonlat0 * std::pow(0.5, i);

            // Make a first order estimate of xy using Jacobian.
            const auto xy = xy0[t] + jac * hLonlat;

            // Get true xy values.
            const auto xyTarget = csProjection.xy(lonlat0[t] + hLonlat);
            double dxy          = Point2::distance(xyTarget, xy);

            // Error should reduce roughly by a factor of four every iteration.
            if (i > 0) {
                const auto tol = 0.02;
                EXPECT(dxy / dxyOld < 0.25 + tol);
            }
            dxyOld = dxy;
        }
    }

    // Create mesh.
    auto mesh = MeshGenerator(meshType).generate(grid);

    // Create functionspace.
    auto functionspace = functionspace::NodeColumns(mesh);

    // Create fields for each component of Jacobian.
    auto fieldSet = FieldSet{};
    fieldSet.add(functionspace.createField<double>(option::name("dx_by_dlambda")));
    fieldSet.add(functionspace.createField<double>(option::name("dx_by_dphi")));
    fieldSet.add(functionspace.createField<double>(option::name("dy_by_dlambda")));
    fieldSet.add(functionspace.createField<double>(option::name("dy_by_dphi")));

    auto dx_dlambda   = array::make_view<double, 1>(fieldSet["dx_by_dlambda"]);
    auto dx_dphi      = array::make_view<double, 1>(fieldSet["dx_by_dphi"]);
    auto dy_dlambda   = array::make_view<double, 1>(fieldSet["dy_by_dlambda"]);
    auto dy_dphi      = array::make_view<double, 1>(fieldSet["dy_by_dphi"]);
    const auto lonlat = array::make_view<double, 2>(functionspace.lonlat());

    for (idx_t i = 0; i < functionspace.size(); ++i) {
        const auto ll  = PointLonLat(lonlat(i, LON), lonlat(i, LAT));
        const auto jac = csProjection.jacobian(ll);
        dx_dlambda(i)  = jac.dx_dlon();
        dx_dphi(i)     = jac.dx_dlat();
        dy_dlambda(i)  = jac.dy_dlon();
        dy_dphi(i)     = jac.dy_dlat();
    }

    // Output mesh.
    const auto gmshConfig =
        util::Config("coordinates", "xy") | util::Config("ghost", true) | util::Config("info", true);
    auto gmsh = output::Gmsh(outputID + "_jacobian.msh", gmshConfig);
    gmsh.write(mesh);
    gmsh.write(fieldSet, functionspace);
}


CASE("cubedsphere_xy_to_alphabeta_test") {
    testProjection("CS-LFR-12", "cubedsphere", "cs_primal");
    testProjection("CS-LFR-12", "cubedsphere_dual", "cs_dual");
}

CASE("cubedsphere_jacobian_test") {
    testJacobian("CS-LFR-24", "cubedsphere", "cs_primal");
    testJacobian("CS-LFR-24", "cubedsphere_dual", "cs_dual");
}

CASE("test_tiles") {
    int resolution(2);
    Grid gEA{"CS-EA-L-" + std::to_string(resolution)};
    Grid gLFR{"CS-LFR-L-" + std::to_string(resolution)};

    using util::Constants;

    util::Config params;
    grid::CubedSphereTiles f("cubedsphere_fv3");
    grid::CubedSphereTiles l("cubedsphere_lfric");

    double cd[2];

    idx_t jn(0);

    std::array<idx_t, 7> EAOffset{0,
                                  resolution * resolution + 1,
                                  2 * resolution * resolution + 2,
                                  3 * resolution * resolution + 2,
                                  4 * resolution * resolution + 2,
                                  5 * resolution * resolution + 2,
                                  6 * resolution * resolution + 2};

    for (auto crd : gEA.lonlat()) {
        atlas::PointLonLat pointLonLat = crd;
        cd[LON]                        = pointLonLat.lon();
        cd[LAT]                        = pointLonLat.lat();

        int t = f.indexFromLonLat(cd);

        gEA.projection().lonlat2xy(crd);
        cd[LON] = crd.lon();
        cd[LAT] = crd.lat();

        int t2 = f.indexFromXY(cd);

        for (std::size_t i = 0; i < 6; ++i) {
            if (jn >= EAOffset[i] && jn < EAOffset[i + 1]) {
                EXPECT(t == static_cast<idx_t>(i));
                EXPECT(t2 == static_cast<idx_t>(i));
            }
        }
        ++jn;
    }

    std::array<idx_t, 7> LFRicOffset{0,
                                     resolution * resolution,
                                     2 * resolution * resolution,
                                     3 * resolution * resolution,
                                     4 * resolution * resolution,
                                     4 * resolution * resolution + (resolution + 1) * (resolution + 1),
                                     6 * resolution * resolution + 2};

    jn = 0;
    for (auto crd : gLFR.lonlat()) {
        atlas::PointLonLat pointLonLat = crd;

        cd[LON] = pointLonLat.lon();
        cd[LAT] = pointLonLat.lat();

        int t3 = l.indexFromLonLat(cd);

        gLFR.projection().lonlat2xy(crd);
        cd[LON] = crd.lon();
        cd[LAT] = crd.lat();

        int t4 = l.indexFromXY(cd);

        for (std::size_t i = 0; i < 6; ++i) {
            if (jn >= LFRicOffset[i] && jn < LFRicOffset[i + 1]) {
                EXPECT(t3 == static_cast<idx_t>(i));
                EXPECT(t4 == static_cast<idx_t>(i));
            }
        }
        ++jn;
    }
}

CASE("test_projection_cubedsphere_xy_latlon") {
    int resolution(12);
    std::vector<std::string> grid_names{"CS-EA-L-" + std::to_string(resolution),
                                        "CS-ED-L-" + std::to_string(resolution)};

    for (std::string& s : grid_names) {
        Grid g{s};
        for (auto crd : g.lonlat()) {
            Point2 lonlat{crd};
            g->projection().lonlat2xy(crd);
            g->projection().xy2lonlat(crd);
            // except for point lonlat (90,82.5) on compiler pgc++
            // we have a maximum error tolerance of 1e-11
            EXPECT_APPROX_EQ(lonlat, crd, 1e-6);
        }
        for (auto crd : g.xy()) {
            Point2 xy{crd};
            g->projection().xy2lonlat(crd);
            g->projection().lonlat2xy(crd);
            EXPECT_APPROX_EQ(xy, crd, 1e-6);
        }
    }
}


}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
