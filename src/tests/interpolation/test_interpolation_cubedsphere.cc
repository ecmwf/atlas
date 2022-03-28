/*
 * (C) Crown Copyright 2022 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */


#include "atlas/field/FieldSet.h"
#include "atlas/functionspace/CellColumns.h"
#include "atlas/functionspace/CubedSphereColumns.h"
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/grid/CubedSphereGrid.h"
#include "atlas/grid/Grid.h"
#include "atlas/grid/Partitioner.h"
#include "atlas/interpolation/Interpolation.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/meshgenerator/MeshGenerator.h"
#include "atlas/output/Gmsh.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/util/Constants.h"
#include "atlas/util/CoordinateEnums.h"
#include "atlas/util/function/VortexRollup.h"

#include "tests/AtlasTestEnvironment.h"


namespace atlas {
namespace test {


// Return (u, v) field with vortex_rollup as the streamfunction.
// This has no physical significance, but it makes a nice swirly field.
std::pair<double, double> vortexField(double lon, double lat) {

    // set hLon and hLat step size.
    const double hLon = 0.0001;
    const double hLat = 0.0001;

    // Get finite differences.

    // Set u.
    const double u =
        (util::function::vortex_rollup(lon, lat + 0.5 * hLat, 1.) -
         util::function::vortex_rollup(lon, lat - 0.5 * hLat, 1.)) / hLat;

    const double v =
        -(util::function::vortex_rollup(lon + 0.5 * hLon, lat, 1.) -
          util::function::vortex_rollup(lon - 0.5 * hLon, lat, 1.)) /
         (hLon * std::cos(lat * util::Constants::degreesToRadians()));

    return std::make_pair(u, v);

}


double dotProd(const Field& a, const Field& b) {
    double prod{};

    const auto aView = array::make_view<double, 1>(a);
    const auto bView = array::make_view<double, 1>(b);


    for (size_t i = 0; i < a.size(); ++i) {
        prod += aView(i) * bView(i);
    }
    mpi::comm().allReduceInPlace(prod, eckit::mpi::Operation::SUM);
    return prod;
}

CASE("cubedsphere_scalar_interpolation") {


    // Create a source cubed sphere grid, mesh and functionspace.
    const auto sourceGrid          = Grid("CS-LFR-24");
    const auto sourceMesh          = MeshGenerator("cubedsphere_dual").generate(sourceGrid);
    const auto sourceFunctionspace = functionspace::NodeColumns(sourceMesh);

    //--------------------------------------------------------------------------
    // Interpolation test.
    //--------------------------------------------------------------------------

    // Populate analytic source field.
    double stDev{};
    auto sourceField = sourceFunctionspace.createField<double>(option::name("test_field"));
    {
        const auto lonlat = array::make_view<double, 2>(sourceFunctionspace.lonlat());
        const auto ghost  = array::make_view<int, 1>(sourceFunctionspace.ghost());
        auto view         = array::make_view<double, 1>(sourceField);
        for (idx_t i = 0; i < sourceFunctionspace.size(); ++i) {
            view(i) = util::function::vortex_rollup(lonlat(i, LON), lonlat(i, LAT), 1.);
            if (!ghost(i)) {
                stDev += view(i) * view(i);
            }
        }
    }
    mpi::comm().allReduceInPlace(stDev, eckit::mpi::Operation::SUM);
    stDev = std::sqrt(stDev / sourceGrid.size());


    // Create target grid, mesh and functionspace.
    const auto partitioner         = grid::MatchingPartitioner(sourceMesh, util::Config("type", "cubedsphere"));
    const auto targetGrid          = Grid("O24");
    const auto targetMesh          = MeshGenerator("structured").generate(targetGrid, partitioner);
    const auto targetFunctionspace = functionspace::NodeColumns(targetMesh);

    // Set up interpolation object.
    const auto scheme = util::Config("type", "cubedsphere-bilinear") | util::Config("adjoint", true);
    const auto interp = Interpolation(scheme, sourceFunctionspace, targetFunctionspace);

    // Interpolate from source to target field.
    auto targetField = targetFunctionspace.createField<double>(option::name("test_field"));
    interp.execute(sourceField, targetField);
    targetField.haloExchange();

    // Make some diagnostic output fields.
    auto errorField = targetFunctionspace.createField<double>(option::name("error_field"));
    auto partField  = targetFunctionspace.createField<int>(option::name("partition"));
    {
        const auto lonlat = array::make_view<double, 2>(targetFunctionspace.lonlat());
        auto targetView   = array::make_view<double, 1>(targetField);
        auto errorView    = array::make_view<double, 1>(errorField);
        auto partView     = array::make_view<int, 1>(partField);
        for (idx_t i = 0; i < targetFunctionspace.size(); ++i) {
            const auto val = util::function::vortex_rollup(lonlat(i, LON), lonlat(i, LAT), 1.);
            errorView(i)   = std::abs((targetView(i) - val) / stDev);
            partView(i)    = mpi::rank();
        }
    }
    partField.haloExchange();

    // Output source mesh.
    const auto gmshConfig =
        util::Config("coordinates", "xyz") | util::Config("ghost", true) | util::Config("info", true);
    const auto sourceGmsh = output::Gmsh("cubedsphere_source.msh", gmshConfig);
    sourceGmsh.write(sourceMesh);
    sourceGmsh.write(FieldSet(sourceField), sourceFunctionspace);

    // Output target mesh.
    const auto targetGmsh = output::Gmsh("cubedsphere_target.msh", gmshConfig);
    targetGmsh.write(targetMesh);
    auto targetFields = FieldSet{};
    targetFields.add(targetField);
    targetFields.add(errorField);
    targetFields.add(partField);
    targetGmsh.write(targetFields, targetFunctionspace);

    //--------------------------------------------------------------------------
    // Adjoint test.
    //--------------------------------------------------------------------------

    // Ensure that the adjoint identity relationship holds.
    auto targetAdjoint = targetFunctionspace.createField<double>(option::name("target adjoint"));
    array::make_view<double, 1>(targetAdjoint).assign(array::make_view<double, 1>(targetField));
    targetAdjoint.adjointHaloExchange();

    auto sourceAdjoint = sourceFunctionspace.createField<double>(option::name("source adjoint"));
    array::make_view<double, 1>(sourceAdjoint).assign(0.);
    interp.execute_adjoint(sourceAdjoint, targetAdjoint);

    const auto yDotY    = dotProd(targetField, targetField);
    const auto xDotXAdj = dotProd(sourceField, sourceAdjoint);

    EXPECT_APPROX_EQ(yDotY / xDotXAdj, 1., 1e-14);
}

CASE("cubedsphere_wind_interpolation") {

    // Create a source cubed sphere grid, mesh and functionspace.
    const auto sourceGrid = CubedSphereGrid("CS-LFR-48");
    const auto sourceMesh = MeshGenerator("cubedsphere_dual").generate(sourceGrid);
    const auto sourceFunctionspace = functionspace::CubedSphereNodeColumns(sourceMesh);

    // Get projection.
    const auto& proj = sourceGrid.cubedSphereProjection();

    // Set wind transform Jacobian.
    const auto windTransform = [&](const PointLonLat& lonlat, idx_t t) {

        // Jacobian = [d(alpha, beta) / d(lon, lat)] * [d(lon, lat) / d(u_hat, v_hat)];
        // Note: u_hat and v_hat are increments along the u and v unit vectors.

        const auto secLat = 1. / std::cos(util::Constants::degreesToRadians() * lonlat.lat());
        return proj.alphabetaJacobian(lonlat, t) * projection::Jacobian{{secLat, 0.}, {0., 1.}};
    };

    // Matrix-vector multiplication helper. y = Ax.
    const auto matMul = [](const projection::Jacobian& A, double x0, double x1){
        const Point2 y = A * Point2{x0, x1};
        return std::make_pair(y[0], y[1]);
    };

    //--------------------------------------------------------------------------
    // Interpolation test.
    //--------------------------------------------------------------------------

    // Populate analytic source field.
    auto sourceFieldSet = FieldSet{};
    sourceFieldSet.add(sourceFunctionspace.createField<double>(option::name("u_orig")));
    sourceFieldSet.add(sourceFunctionspace.createField<double>(option::name("v_orig")));
    sourceFieldSet.add(sourceFunctionspace.createField<double>(option::name("v_alpha")));
    sourceFieldSet.add(sourceFunctionspace.createField<double>(option::name("v_beta")));
    {
        const auto lonlat = array::make_view<double, 2>(sourceFunctionspace.lonlat());
        auto u = array::make_view<double, 1>(sourceFieldSet["u_orig"]);
        auto v = array::make_view<double, 1>(sourceFieldSet["v_orig"]);
        auto vAlpha = array::make_view<double, 1>(sourceFieldSet["v_alpha"]);
        auto vBeta = array::make_view<double, 1>(sourceFieldSet["v_beta"]);


        // In practice, the (u, v) field needs to be halo exchanged *before* the
        // wind transform. Then the transform is applied to the entire field,
        // *including* the halo.

        sourceFunctionspace.parallel_for(util::Config("include_halo", true),
            [&](idx_t idx, idx_t t, idx_t i, idx_t j){

                // Get lonlat
                const auto ll = PointLonLat(lonlat(idx, LON), lonlat(idx, LAT));

                // Set (u, v) wind
                std::tie(u(idx), v(idx)) = vortexField(ll.lon(), ll.lat());

                // Get wind transform jacobian.
                auto jac = windTransform(ll, t);

                // Transform wind.
                std::tie(vAlpha(idx), vBeta(idx)) = matMul(jac, u(idx), v(idx));

            });
    }

    // Create target grid, mesh and functionspace.
    const auto partitioner = grid::MatchingPartitioner(sourceMesh, util::Config("type", "cubedsphere"));
    const auto targetGrid = Grid("O48");
    const auto targetMesh = MeshGenerator("structured").generate(targetGrid, partitioner);
    const auto targetFunctionspace = functionspace::NodeColumns(targetMesh);

    // Set up interpolation object.
    // Note: We have to disable the source field halo exhange in the
    // interpolation execute and excute_adjoint methods. If left on, the halo
    // exchange will corrupt the transformed wind field.
    const auto scheme = util::Config("type", "cubedsphere-bilinear") |
                        util::Config("adjoint", true) |
                        util::Config("halo_exchange", false);
    const auto interp = Interpolation(scheme, sourceFunctionspace, targetFunctionspace);

    // Make target fields.
    auto targetFieldSet = FieldSet{};
    targetFieldSet.add(targetFunctionspace.createField<double>(option::name("u_orig")));
    targetFieldSet.add(targetFunctionspace.createField<double>(option::name("v_orig")));
    targetFieldSet.add(targetFunctionspace.createField<double>(option::name("v_alpha")));
    targetFieldSet.add(targetFunctionspace.createField<double>(option::name("v_beta")));

    // Interpolate from source to target fields.
    array::make_view<double, 1>(targetFieldSet["v_alpha"]).assign(0.);
    array::make_view<double, 1>(targetFieldSet["v_beta"]).assign(0.);
    interp.execute(sourceFieldSet, targetFieldSet);

    // Make new (u, v) fields from (v_alpha, v_beta)
    targetFieldSet.add(targetFunctionspace.createField<double>(option::name("error_field_0")));
    targetFieldSet.add(targetFunctionspace.createField<double>(option::name("error_field_1")));
    targetFieldSet.add(targetFunctionspace.createField<double>(option::name("u_new")));
    targetFieldSet.add(targetFunctionspace.createField<double>(option::name("v_new")));
    {
        const auto lonlat = array::make_view<double, 2>(targetFunctionspace.lonlat());
        const auto ghost = array::make_view<int, 1>(targetFunctionspace.ghost());
        auto u = array::make_view<double, 1>(targetFieldSet["u_new"]);
        auto v = array::make_view<double, 1>(targetFieldSet["v_new"]);
        const auto uOrig = array::make_view<double, 1>(targetFieldSet["u_orig"]);
        const auto vOrig = array::make_view<double, 1>(targetFieldSet["v_orig"]);
        const auto vAlpha = array::make_view<double, 1>(targetFieldSet["v_alpha"]);
        const auto vBeta = array::make_view<double, 1>(targetFieldSet["v_beta"]);
        auto error0 = array::make_view<double, 1>(targetFieldSet["error_field_0"]);
        auto error1 = array::make_view<double, 1>(targetFieldSet["error_field_1"]);
        for (idx_t idx = 0; idx < targetFunctionspace.size(); ++idx) {

            if (!ghost(idx)) {
                const auto ll = PointLonLat(lonlat(idx, LON), lonlat(idx, LAT));
                const idx_t t = proj.getCubedSphereTiles().indexFromLonLat(ll);

                // Get inverse wind transform jacobian.
                auto jac = windTransform(ll, t).inverse();

                // Transform wind.
                std::tie(u(idx), v(idx)) = matMul(jac, vAlpha(idx), vBeta(idx));

                // Get error.
                const auto uvTarget = vortexField(ll.lon(), ll.lat());

                error0(idx) = Point2::distance(Point2(uvTarget.first, uvTarget.second), Point2(uOrig(idx), vOrig(idx)));
                error1(idx) = Point2::distance(Point2(uvTarget.first, uvTarget.second), Point2(u(idx), v(idx)));
            }
        }
    }
    targetFieldSet.haloExchange();

    // Output source mesh.
    const auto gmshConfig =
        util::Config("coordinates", "xyz") | util::Config("ghost", true) | util::Config("info", true);
    const auto sourceGmsh = output::Gmsh("cubedsphere_vec_source.msh", gmshConfig);
    sourceGmsh.write(sourceMesh);
    sourceGmsh.write(sourceFieldSet, sourceFunctionspace);

    // Output target mesh.
    const auto targetGmsh = output::Gmsh("cubedsphere_vec_target.msh", gmshConfig);
    targetGmsh.write(targetMesh);
    targetGmsh.write(targetFieldSet, targetFunctionspace);


    //--------------------------------------------------------------------------
    // Adjoint test.
    //--------------------------------------------------------------------------

    // Ensure that the adjoint identity relationship holds.
    targetFieldSet.add(targetFunctionspace.createField<double>(option::name("u_adjoint")));
    targetFieldSet.add(targetFunctionspace.createField<double>(option::name("v_adjoint")));
    array::make_view<double, 1>(targetFieldSet["u_adjoint"]).assign(array::make_view<double, 1>(targetFieldSet["u_new"]));
    array::make_view<double, 1>(targetFieldSet["v_adjoint"]).assign(array::make_view<double, 1>(targetFieldSet["v_new"]));

    // Adjoint of target halo exhange.
    targetFieldSet["u_adjoint"].adjointHaloExchange();
    targetFieldSet["v_adjoint"].adjointHaloExchange();

    // Adjoint of inverse wind transform.
    targetFieldSet.add(targetFunctionspace.createField<double>(option::name("v_alpha_adjoint")));
    targetFieldSet.add(targetFunctionspace.createField<double>(option::name("v_beta_adjoint")));
    {
        const auto lonlat = array::make_view<double, 2>(targetFunctionspace.lonlat());
        const auto ghost = array::make_view<int, 1>(targetFunctionspace.ghost());
        const auto uAdj = array::make_view<double, 1>(targetFieldSet["u_adjoint"]);
        const auto vAdj = array::make_view<double, 1>(targetFieldSet["v_adjoint"]);
        auto vAlphaAdj = array::make_view<double, 1>(targetFieldSet["v_alpha_adjoint"]);
        auto vBetaAdj = array::make_view<double, 1>(targetFieldSet["v_beta_adjoint"]);
        vAlphaAdj.assign(0.);
        vBetaAdj.assign(0.);
        for (idx_t idx = 0; idx < targetFunctionspace.size(); ++idx) {

            if (!ghost(idx)) {
                const auto ll = PointLonLat(lonlat(idx, LON), lonlat(idx, LAT));
                const idx_t t = proj.getCubedSphereTiles().indexFromLonLat(ll);

                // Get adjoint of inverse wind transform jacobian.
                auto jac = windTransform(ll, t).inverse().transpose();

                // Transform wind.
                std::tie(vAlphaAdj(idx), vBetaAdj(idx)) = matMul(jac, uAdj(idx), vAdj(idx));
            }
        }
    }

    // Adjoint of interpolation.
    sourceFieldSet.add(sourceFunctionspace.createField<double>(option::name("v_alpha_adjoint")));
    sourceFieldSet.add(sourceFunctionspace.createField<double>(option::name("v_beta_adjoint")));
    array::make_view<double, 1>(sourceFieldSet["v_alpha_adjoint"]).assign(0.);
    array::make_view<double, 1>(sourceFieldSet["v_beta_adjoint"]).assign(0.);
    interp.execute_adjoint(sourceFieldSet["v_alpha_adjoint"], targetFieldSet["v_alpha_adjoint"]);
    interp.execute_adjoint(sourceFieldSet["v_beta_adjoint"], targetFieldSet["v_beta_adjoint"]);

    // Adjoint of wind transform.
    sourceFieldSet.add(sourceFunctionspace.createField<double>(option::name("u_adjoint")));
    sourceFieldSet.add(sourceFunctionspace.createField<double>(option::name("v_adjoint")));
    {
        const auto lonlat = array::make_view<double, 2>(sourceFunctionspace.lonlat());
        auto uAdj = array::make_view<double, 1>(sourceFieldSet["u_adjoint"]);
        auto vAdj = array::make_view<double, 1>(sourceFieldSet["v_adjoint"]);
        uAdj.assign(0.);
        vAdj.assign(0.);
        const auto vAlphaAdj = array::make_view<double, 1>(sourceFieldSet["v_alpha_adjoint"]);
        const auto vBetaAdj = array::make_view<double, 1>(sourceFieldSet["v_beta_adjoint"]);

        sourceFunctionspace.parallel_for(util::Config("include_halo", true),
            [&](idx_t idx, idx_t t, idx_t i, idx_t j){

                // Get lonlat
                const auto ll = PointLonLat(lonlat(idx, LON), lonlat(idx, LAT));

                // Get adjoint of wind transform jacobian.
                auto jac = windTransform(ll, t).transpose();

                // Transform wind.
                std::tie(uAdj(idx), vAdj(idx)) = matMul(jac, vAlphaAdj(idx), vBetaAdj(idx));
            });
    }

    // Check dot products.
    const auto yDotY = dotProd(targetFieldSet["u_new"], targetFieldSet["u_new"]) +
                       dotProd(targetFieldSet["v_new"], targetFieldSet["v_new"]);

    const auto xDotXAdj = dotProd(sourceFieldSet["u_orig"], sourceFieldSet["u_adjoint"]) +
                          dotProd(sourceFieldSet["v_orig"], sourceFieldSet["v_adjoint"]);

    EXPECT_APPROX_EQ(yDotY / xDotXAdj, 1., 1e-14);


}


}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
