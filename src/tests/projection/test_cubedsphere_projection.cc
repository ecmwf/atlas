/*
 * (C) Crown Copyright 2022 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "atlas/array/MakeView.h"
#include "atlas/grid.h"
#include "atlas/grid/CubedSphereGrid.h"
#include "atlas/mesh.h"
#include "atlas/meshgenerator.h"
#include "atlas/projection/detail/CubedSphereProjectionBase.h"
#include "atlas/output/Gmsh.h"
#include "atlas/util/Matrix.h"
#include "tests/AtlasTestEnvironment.h"
namespace atlas {
namespace test {

// Small number relative to 360.
constexpr double epsilon = std::numeric_limits<double>::epsilon() * 360.;

void testProjection(const std::string& gridType, const std::string& meshType,
                    const std::string& outputID) {

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
                    << csProjection.getCubedSphereTiles().tileJacobian(i) << std::endl << std::endl;
    }

    // Create mesh.
    auto mesh = MeshGenerator(meshType, util::Config("halo", 3)).generate(grid);

    // Output mesh.
    const auto gmshConfig =
        util::Config("coordinates", "xy") | util::Config("ghost", true) | util::Config("info", true);
    auto gmsh = output::Gmsh(outputID + "_before.msh", gmshConfig);
    gmsh.write(mesh);



    // "Hack" mesh xy coordinates.
    auto xyView = array::make_view<double, 2>(mesh.nodes().xy());
    const auto tijView = array::make_view<idx_t, 2>(mesh.nodes().field("tij"));
    for (idx_t i = 0; i < mesh.nodes().size(); ++i) {
        const auto t = tijView(i, 0);

        const auto xy = PointXY(xyView(i, XX), xyView(i, YY));
        const auto alphabeta = csProjection.xy2alphabeta(xy, t);

        // Inverse function is degenerate when abs(alpha) == abs(beta) and
        // abs(alpha) > 45.
        const bool degenerate = std::abs(alphabeta[0]) > 45. &&
            approx_eq(std::abs(alphabeta[0]), std::abs(alphabeta[1]));

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



CASE("cubedsphere_xy_to_alphabeta_test") {

    testProjection("CS-LFR-13", "cubedsphere", "cs_primal");
    testProjection("CS-LFR-13", "cubedsphere_dual", "cs_dual");

}


}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
