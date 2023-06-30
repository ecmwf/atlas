/*
 * (C) Crown Copyright 2022 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */


#include "atlas/array/helpers/ArrayForEach.h"
#include "atlas/field/FieldSet.h"
#include "atlas/functionspace/CellColumns.h"
#include "atlas/functionspace/CubedSphereColumns.h"
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/functionspace/PointCloud.h"
#include "atlas/functionspace/StructuredColumns.h"
#include "atlas/grid/CubedSphereGrid.h"
#include "atlas/grid/Distribution.h"
#include "atlas/grid/Grid.h"
#include "atlas/grid/Iterator.h"
#include "atlas/grid/Partitioner.h"
#include "atlas/interpolation/Interpolation.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/meshgenerator/MeshGenerator.h"
#include "atlas/output/Gmsh.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/redistribution/Redistribution.h"
#include "atlas/util/Constants.h"
#include "atlas/util/CoordinateEnums.h"
#include "atlas/util/function/VortexRollup.h"

#include "tests/AtlasTestEnvironment.h"


namespace atlas {
namespace test {

struct CubedSphereInterpolationFixture {
    atlas::Grid sourceGrid_ = Grid("CS-LFR-24");
    atlas::MeshGenerator sourceMeshGenerator_ = MeshGenerator("cubedsphere_dual");
    atlas::grid::Partitioner sourcePartitioner_ = sourceGrid_.partitioner();
    atlas::Mesh sourceMesh_ = sourceMeshGenerator_.generate(sourceGrid_, sourcePartitioner_);
    atlas::FunctionSpace sourceFunctionSpace_ = functionspace::NodeColumns(sourceMesh_);
    atlas::Grid targetGrid_ = Grid("O24");
    atlas::MeshGenerator targetMeshGenerator_ =  MeshGenerator("structured");
    atlas::grid::Partitioner targetPartitioner_ =
        grid::MatchingPartitioner(sourceMesh_, util::Config("type", "cubedsphere"));
    atlas::Mesh targetMesh_ = targetMeshGenerator_.generate(targetGrid_, targetPartitioner_);
    atlas::FunctionSpace targetFunctionSpace_ = functionspace::NodeColumns(targetMesh_);
};

void gmshOutput(const std::string& fileName, const FieldSet& fieldSet) {


    const auto& functionSpace = fieldSet[0].functionspace();
    const auto& mesh = functionspace::NodeColumns(functionSpace).mesh();

    const auto gmshConfig =
        util::Config("coordinates", "xyz") | util::Config("ghost", true) | util::Config("info", true);
    const auto gmsh = output::Gmsh(fileName, gmshConfig);
    gmsh.write(mesh);
    gmsh.write(fieldSet, functionSpace);
}

atlas::util::Config interpScheme() {
    return util::Config("type", "cubedsphere-bilinear") | util::Config("adjoint", true) |
           util::Config("halo_exchange", false);
}

atlas::Field generateScalarField(const FunctionSpace& functionspace) {

    auto field = functionspace.createField<double>(option::name("test_field"));
    auto fieldView = array::make_view<double, 1>(field);
    const auto lonlLatView = array::make_view<double, 2>(functionspace.lonlat());

    using for_each = array::helpers::ArrayForEach<0>;

    for_each::apply(std::make_tuple(fieldView, lonlLatView), [&](auto&& fieldElem, auto&& lonLat){
        fieldElem = util::function::vortex_rollup(lonLat(LON), lonLat(LAT), 1.);
    });
    return field;
};


// Return (u, v) field with vortex_rollup as the streamfunction.
// This has no physical significance, but it makes a nice swirly field.
std::pair<double, double> vortexVectorField(double lon, double lat) {

    // set hLon and hLat step size.
    const double hLon = 0.0001;
    const double hLat = 0.0001;

    // Get finite differences.

    // Set u.
    const double u = (util::function::vortex_rollup(lon, lat + 0.5 * hLat, 1.) -
                      util::function::vortex_rollup(lon, lat - 0.5 * hLat, 1.)) /
                     hLat;

    const double v = -(util::function::vortex_rollup(lon + 0.5 * hLon, lat, 1.) -
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

    const auto fixture = CubedSphereInterpolationFixture{};

    //--------------------------------------------------------------------------
    // Interpolation test.
    //--------------------------------------------------------------------------

    // Populate analytic source field.
    double stDev{};
    auto sourceField = fixture.sourceFunctionSpace_.createField<double>(option::name("test_field"));
    {
        const auto lonlat = array::make_view<double, 2>(fixture.sourceFunctionSpace_.lonlat());
        const auto ghost  = array::make_view<int, 1>(fixture.sourceFunctionSpace_.ghost());
        auto view         = array::make_view<double, 1>(sourceField);
        for (idx_t i = 0; i < fixture.sourceFunctionSpace_.size(); ++i) {
            view(i) = util::function::vortex_rollup(lonlat(i, LON), lonlat(i, LAT), 1.);
            if (!ghost(i)) {
                stDev += view(i) * view(i);
            }
        }
    }
    mpi::comm().allReduceInPlace(stDev, eckit::mpi::Operation::SUM);
    stDev = std::sqrt(stDev / fixture.sourceGrid_.size());

    // Set up interpolation object.
    const auto interp = Interpolation(interpScheme(), fixture.sourceFunctionSpace_, fixture.targetFunctionSpace_);

    // Interpolate from source to target field.
    auto targetField = fixture.targetFunctionSpace_.createField<double>(option::name("test_field"));
    interp.execute(sourceField, targetField);
    targetField.haloExchange();

    // Make some diagnostic output fields.
    auto errorField = fixture.targetFunctionSpace_.createField<double>(option::name("error_field"));
    auto partField  = fixture.targetFunctionSpace_.createField<int>(option::name("partition"));
    {
        const auto lonlat = array::make_view<double, 2>(fixture.targetFunctionSpace_.lonlat());
        auto targetView   = array::make_view<double, 1>(targetField);
        auto errorView    = array::make_view<double, 1>(errorField);
        auto partView     = array::make_view<int, 1>(partField);
        for (idx_t i = 0; i < fixture.targetFunctionSpace_.size(); ++i) {
            const auto val = util::function::vortex_rollup(lonlat(i, LON), lonlat(i, LAT), 1.);
            errorView(i)   = std::abs((targetView(i) - val) / stDev);
            partView(i)    = mpi::rank();
        }
    }
    partField.haloExchange();

    gmshOutput("cubedsphere_source.msh", FieldSet(sourceField));

    auto targetFields = FieldSet{};
    targetFields.add(targetField);
    targetFields.add(errorField);
    targetFields.add(partField);
    gmshOutput("cubedsphere_target.msh", targetFields);

    //--------------------------------------------------------------------------
    // Adjoint test.
    //--------------------------------------------------------------------------

    // Ensure that the adjoint identity relationship holds.
    auto targetAdjoint = fixture.targetFunctionSpace_.createField<double>(option::name("target adjoint"));
    array::make_view<double, 1>(targetAdjoint).assign(array::make_view<double, 1>(targetField));
    targetAdjoint.adjointHaloExchange();

    auto sourceAdjoint = fixture.sourceFunctionSpace_.createField<double>(option::name("source adjoint"));
    array::make_view<double, 1>(sourceAdjoint).assign(0.);
    interp.execute_adjoint(sourceAdjoint, targetAdjoint);

    const auto yDotY    = dotProd(targetField, targetField);
    const auto xDotXAdj = dotProd(sourceField, sourceAdjoint);

    EXPECT_APPROX_EQ(yDotY / xDotXAdj, 1., 1e-14);
}

CASE("cubedsphere_wind_interpolation") {

    const auto fixture = CubedSphereInterpolationFixture{};

    // Get projection.
    const auto& proj = CubedSphereGrid(fixture.sourceGrid_).cubedSphereProjection();

    // Set wind transform Jacobian.
    const auto windTransform = [&](const PointLonLat& lonlat, idx_t t) {
        // Jacobian = [d(alpha, beta) / d(lon, lat)] * [d(lon, lat) / d(u_hat, v_hat)];
        // Note: u_hat and v_hat are increments along the u and v unit vectors.

        const auto secLat = 1. / std::cos(util::Constants::degreesToRadians() * lonlat.lat());
        return proj.alphabetaJacobian(lonlat, t) * projection::Jacobian{{secLat, 0.}, {0., 1.}};
    };

    // Matrix-vector multiplication helper. y = Ax.
    const auto matMul = [](const projection::Jacobian& A, double x0, double x1) {
        const Point2 y = A * Point2{x0, x1};
        return std::make_pair(y[0], y[1]);
    };

    //--------------------------------------------------------------------------
    // Interpolation test.
    //--------------------------------------------------------------------------

    // Populate analytic source field.
    auto sourceFieldSet = FieldSet{};
    sourceFieldSet.add(fixture.sourceFunctionSpace_.createField<double>(option::name("u_orig")));
    sourceFieldSet.add(fixture.sourceFunctionSpace_.createField<double>(option::name("v_orig")));
    sourceFieldSet.add(fixture.sourceFunctionSpace_.createField<double>(option::name("v_alpha")));
    sourceFieldSet.add(fixture.sourceFunctionSpace_.createField<double>(option::name("v_beta")));
    {
        const auto lonlat = array::make_view<double, 2>(fixture.sourceFunctionSpace_.lonlat());
        auto u            = array::make_view<double, 1>(sourceFieldSet["u_orig"]);
        auto v            = array::make_view<double, 1>(sourceFieldSet["v_orig"]);
        auto vAlpha       = array::make_view<double, 1>(sourceFieldSet["v_alpha"]);
        auto vBeta        = array::make_view<double, 1>(sourceFieldSet["v_beta"]);


        // In practice, the (u, v) field needs to be halo exchanged *before* the
        // wind transform. Then the transform is applied to the entire field,
        // *including* the halo.

        functionspace::CubedSphereNodeColumns(fixture.sourceFunctionSpace_).parallel_for(
            util::Config("include_halo", true), [&](idx_t idx, idx_t t, idx_t i, idx_t j) {
            // Get lonlat
            const auto ll = PointLonLat(lonlat(idx, LON), lonlat(idx, LAT));

            // Set (u, v) wind
            std::tie(u(idx), v(idx)) = vortexVectorField(ll.lon(), ll.lat());

            // Get wind transform jacobian.
            auto jac = windTransform(ll, t);

            // Transform wind.
            std::tie(vAlpha(idx), vBeta(idx)) = matMul(jac, u(idx), v(idx));
        });
    }

    // Set up interpolation object.
    // Note: We have to disable the source field halo exhange in the
    // interpolation execute and excute_adjoint methods. If left on, the halo
    // exchange will corrupt the transformed wind field.
    const auto interp = Interpolation(interpScheme(), fixture.sourceFunctionSpace_, fixture.targetFunctionSpace_);

    // Make target fields.
    auto targetFieldSet = FieldSet{};
    targetFieldSet.add(fixture.targetFunctionSpace_.createField<double>(option::name("u_orig")));
    targetFieldSet.add(fixture.targetFunctionSpace_.createField<double>(option::name("v_orig")));
    targetFieldSet.add(fixture.targetFunctionSpace_.createField<double>(option::name("v_alpha")));
    targetFieldSet.add(fixture.targetFunctionSpace_.createField<double>(option::name("v_beta")));

    // Interpolate from source to target fields.
    array::make_view<double, 1>(targetFieldSet["v_alpha"]).assign(0.);
    array::make_view<double, 1>(targetFieldSet["v_beta"]).assign(0.);
    interp.execute(sourceFieldSet, targetFieldSet);

    // Make new (u, v) fields from (v_alpha, v_beta)
    targetFieldSet.add(fixture.targetFunctionSpace_.createField<double>(option::name("error_field_0")));
    targetFieldSet.add(fixture.targetFunctionSpace_.createField<double>(option::name("error_field_1")));
    targetFieldSet.add(fixture.targetFunctionSpace_.createField<double>(option::name("u_new")));
    targetFieldSet.add(fixture.targetFunctionSpace_.createField<double>(option::name("v_new")));
    {
        const auto lonlat = array::make_view<double, 2>(fixture.targetFunctionSpace_.lonlat());
        const auto ghost  = array::make_view<int, 1>(fixture.targetFunctionSpace_.ghost());
        auto u            = array::make_view<double, 1>(targetFieldSet["u_new"]);
        auto v            = array::make_view<double, 1>(targetFieldSet["v_new"]);
        const auto uOrig  = array::make_view<double, 1>(targetFieldSet["u_orig"]);
        const auto vOrig  = array::make_view<double, 1>(targetFieldSet["v_orig"]);
        const auto vAlpha = array::make_view<double, 1>(targetFieldSet["v_alpha"]);
        const auto vBeta  = array::make_view<double, 1>(targetFieldSet["v_beta"]);
        auto error0       = array::make_view<double, 1>(targetFieldSet["error_field_0"]);
        auto error1       = array::make_view<double, 1>(targetFieldSet["error_field_1"]);
        const auto& tVec  = interp.target()->metadata().getIntVector("tile index");

        for (idx_t idx = 0; idx < fixture.targetFunctionSpace_.size(); ++idx) {
            if (!ghost(idx)) {
                const auto ll = PointLonLat(lonlat(idx, LON), lonlat(idx, LAT));
                const idx_t t = tVec[idx];

                // Get inverse wind transform jacobian.
                auto jac = windTransform(ll, t).inverse();

                // Transform wind.
                std::tie(u(idx), v(idx)) = matMul(jac, vAlpha(idx), vBeta(idx));

                // Get error.
                const auto uvTarget = vortexVectorField(ll.lon(), ll.lat());

                error0(idx) = Point2::distance(Point2(uvTarget.first, uvTarget.second), Point2(uOrig(idx), vOrig(idx)));
                error1(idx) = Point2::distance(Point2(uvTarget.first, uvTarget.second), Point2(u(idx), v(idx)));
            }
        }
    }
    targetFieldSet.haloExchange();

    gmshOutput("cubedsphere_vec_source.msh", sourceFieldSet);
    gmshOutput("cubedsphere_vec_target.msh", targetFieldSet);

    //--------------------------------------------------------------------------
    // Adjoint test.
    //--------------------------------------------------------------------------

    // Ensure that the adjoint identity relationship holds.
    targetFieldSet.add(fixture.targetFunctionSpace_.createField<double>(option::name("u_adjoint")));
    targetFieldSet.add(fixture.targetFunctionSpace_.createField<double>(option::name("v_adjoint")));
    array::make_view<double, 1>(targetFieldSet["u_adjoint"])
        .assign(array::make_view<double, 1>(targetFieldSet["u_new"]));
    array::make_view<double, 1>(targetFieldSet["v_adjoint"])
        .assign(array::make_view<double, 1>(targetFieldSet["v_new"]));

    // Adjoint of target halo exhange.
    targetFieldSet["u_adjoint"].adjointHaloExchange();
    targetFieldSet["v_adjoint"].adjointHaloExchange();

    // Adjoint of inverse wind transform.
    targetFieldSet.add(fixture.targetFunctionSpace_.createField<double>(option::name("v_alpha_adjoint")));
    targetFieldSet.add(fixture.targetFunctionSpace_.createField<double>(option::name("v_beta_adjoint")));
    {
        const auto lonlat = array::make_view<double, 2>(fixture.targetFunctionSpace_.lonlat());
        const auto ghost  = array::make_view<int, 1>(fixture.targetFunctionSpace_.ghost());
        const auto uAdj   = array::make_view<double, 1>(targetFieldSet["u_adjoint"]);
        const auto vAdj   = array::make_view<double, 1>(targetFieldSet["v_adjoint"]);
        auto vAlphaAdj    = array::make_view<double, 1>(targetFieldSet["v_alpha_adjoint"]);
        auto vBetaAdj     = array::make_view<double, 1>(targetFieldSet["v_beta_adjoint"]);
        const auto& tVec  = interp.target()->metadata().getIntVector("tile index");
        vAlphaAdj.assign(0.);
        vBetaAdj.assign(0.);
        for (idx_t idx = 0; idx < fixture.targetFunctionSpace_.size(); ++idx) {
            if (!ghost(idx)) {
                const auto ll = PointLonLat(lonlat(idx, LON), lonlat(idx, LAT));
                const idx_t t = tVec[idx];

                // Get adjoint of inverse wind transform jacobian.
                auto jac = windTransform(ll, t).inverse().transpose();

                // Transform wind.
                std::tie(vAlphaAdj(idx), vBetaAdj(idx)) = matMul(jac, uAdj(idx), vAdj(idx));
            }
        }
    }

    // Adjoint of interpolation.
    sourceFieldSet.add(fixture.sourceFunctionSpace_.createField<double>(option::name("v_alpha_adjoint")));
    sourceFieldSet.add(fixture.sourceFunctionSpace_.createField<double>(option::name("v_beta_adjoint")));
    array::make_view<double, 1>(sourceFieldSet["v_alpha_adjoint"]).assign(0.);
    array::make_view<double, 1>(sourceFieldSet["v_beta_adjoint"]).assign(0.);
    interp.execute_adjoint(sourceFieldSet["v_alpha_adjoint"], targetFieldSet["v_alpha_adjoint"]);
    interp.execute_adjoint(sourceFieldSet["v_beta_adjoint"], targetFieldSet["v_beta_adjoint"]);

    // Adjoint of wind transform.
    sourceFieldSet.add(fixture.sourceFunctionSpace_.createField<double>(option::name("u_adjoint")));
    sourceFieldSet.add(fixture.sourceFunctionSpace_.createField<double>(option::name("v_adjoint")));
    {
        const auto lonlat = array::make_view<double, 2>(fixture.sourceFunctionSpace_.lonlat());
        auto uAdj         = array::make_view<double, 1>(sourceFieldSet["u_adjoint"]);
        auto vAdj         = array::make_view<double, 1>(sourceFieldSet["v_adjoint"]);
        uAdj.assign(0.);
        vAdj.assign(0.);
        const auto vAlphaAdj = array::make_view<double, 1>(sourceFieldSet["v_alpha_adjoint"]);
        const auto vBetaAdj  = array::make_view<double, 1>(sourceFieldSet["v_beta_adjoint"]);

        functionspace::CubedSphereNodeColumns(fixture.sourceFunctionSpace_).parallel_for(
            util::Config("include_halo", true), [&](idx_t idx, idx_t t, idx_t i, idx_t j) {
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

CASE("cubedsphere_node_columns_to_structured_columns") {

    const auto fixture = CubedSphereInterpolationFixture{};

    // Can't (easily) redistribute directly from a cubedsphere functionspace to a structured columns.
    // Solution is to build two intermediate PointClouds functions spaces, and copy fields in and out.


    const auto targetStructuredColumns = functionspace::StructuredColumns(fixture.targetGrid_, grid::Partitioner("equal_regions"));
    const auto& targetCubedSphereParitioner = fixture.targetPartitioner_;
    const auto targetNativePartitioner = grid::Partitioner(targetStructuredColumns.distribution());

    // This should be a PointCloud constructor.
    const auto makePointCloud = [](const Grid& grid, const grid::Partitioner partitioner) {

        const auto distribution = grid::Distribution(grid, partitioner);

        auto lonLats = std::vector<PointXY>{};
        auto idx = gidx_t{0};
        for (const auto& lonLat : grid.lonlat()) {
            if (distribution.partition(idx++) == mpi::rank()) {
                lonLats.emplace_back(lonLat.data());
            }
        }
        return functionspace::PointCloud(lonLats);
    };

    const auto targetNativePointCloud = makePointCloud(fixture.targetGrid_, targetNativePartitioner);
    const auto targetCubedSpherePointCloud = makePointCloud(fixture.targetGrid_, targetCubedSphereParitioner);


    // Populate analytic source field.
    auto sourceField = generateScalarField(fixture.sourceFunctionSpace_);

    // Interpolate from source field to targetCubedSpherePointCloud field.
    const auto scheme = util::Config("type", "cubedsphere-bilinear") | util::Config("adjoint", true) |
                        util::Config("halo_exchange", false);
    const auto interp = Interpolation(scheme, fixture.sourceFunctionSpace_, targetCubedSpherePointCloud);
    auto targetCubedSphereField = targetCubedSpherePointCloud.createField<double>(option::name("test_field"));
    interp.execute(sourceField, targetCubedSphereField);

    // Redistribute from targetCubedSpherePointCloud to targetNativePointCloud
    const auto redist = Redistribution(targetCubedSpherePointCloud, targetNativePointCloud);
    auto targetNativeField = targetNativePointCloud.createField<double>(option::name("test_field"));
    redist.execute(targetCubedSphereField, targetNativeField);

    // copy temp field to target field.
    auto targetField = targetStructuredColumns.createField<double>(option::name("test_field"));
    array::make_view<double, 1>(targetField).assign(array::make_view<double, 1>(targetNativeField));

    // Done. Tidy up and write output.

    targetField.haloExchange();

    gmshOutput("cubedsphere_to_structured_cols_source.msh", FieldSet{sourceField});

    // gmsh needs a target mesh...
    const auto mesh = MeshGenerator("structured").generate(fixture.targetGrid_, targetNativePartitioner);
    const auto gmshConfig =
        util::Config("coordinates", "xyz") | util::Config("ghost", true) | util::Config("info", true);
    const auto targetGmsh = output::Gmsh("cubedsphere_to_structured_cols_target.msh", gmshConfig);
    targetGmsh.write(mesh);
    targetGmsh.write(targetField, functionspace::NodeColumns(mesh));
}

CASE("cubedsphere_paritioning_consistency") {

    const auto fixture = CubedSphereInterpolationFixture{};

    const auto makeFunctionSpace = [](const Grid& grid,
                                      const MeshGenerator& meshGenerator,
                                      const grid::Partitioner& partitioner) -> FunctionSpace {
        const auto mesh = meshGenerator.generate(grid, partitioner);
        return functionspace::NodeColumns(mesh);
    };

    const auto& sourceFunctionSpaceCubedSphere = fixture.sourceFunctionSpace_;
    const auto sourceFunctionSpaceEqualRegions = makeFunctionSpace(fixture.sourceGrid_,
                                                                   fixture.sourceMeshGenerator_,
                                                                   grid::Partitioner("equal_regions"));

    const auto& targetFunctionSpaceCubedSphere = fixture.targetFunctionSpace_;

    const auto& targetMeshEqualRegions = functionspace::NodeColumns(sourceFunctionSpaceEqualRegions).mesh();
    const auto targetFunctionSpaceEqualRegions = makeFunctionSpace(fixture.targetGrid_,
                                                                   fixture.targetMeshGenerator_,
                                                                   grid::MatchingPartitioner(
                                                                   targetMeshEqualRegions, util::Config("type", "cubedsphere")));

    const auto interpCubedSphere = Interpolation(interpScheme(), sourceFunctionSpaceCubedSphere, targetFunctionSpaceCubedSphere);
    const auto interpEqualRegions = Interpolation(interpScheme(), sourceFunctionSpaceEqualRegions, targetFunctionSpaceEqualRegions);

    const auto sourceFieldCubedSphere = generateScalarField(sourceFunctionSpaceCubedSphere);
    const auto sourceFieldEqualRegions = generateScalarField(sourceFunctionSpaceEqualRegions);
    auto targetFieldCubedSphere = targetFunctionSpaceCubedSphere.createField<double>(option::name("test field"));
    auto targetFieldEqualRegions = targetFunctionSpaceEqualRegions.createField<double>(option::name("test field"));

    interpCubedSphere.execute(sourceFieldCubedSphere, targetFieldCubedSphere);
    interpEqualRegions.execute(sourceFieldEqualRegions, targetFieldEqualRegions);

    auto targetFieldCubedSphereGlobal = targetFunctionSpaceCubedSphere.createField<double>(option::name("test field global") | option::global(0));
    auto targetFieldEqualRegionsGlobal = targetFunctionSpaceEqualRegions.createField<double>(option::name("test field global") | option::global(0));

    targetFunctionSpaceCubedSphere.gather(targetFieldCubedSphere, targetFieldCubedSphereGlobal);
    targetFunctionSpaceEqualRegions.gather(targetFieldEqualRegions, targetFieldEqualRegionsGlobal);



    auto targetFieldDiffGlobal = targetFunctionSpaceCubedSphere.createField<double>(option::name("test field diff") | option::global(0));
    auto targetFieldDiff = targetFunctionSpaceCubedSphere.createField<double>(option::name("test field diff"));
    auto maxDiff = 0.;

    const auto targetViewCubedSphere = array::make_view<double, 1>(targetFieldCubedSphereGlobal);
    const auto targetViewEqualRegions = array::make_view<double, 1>(targetFieldEqualRegionsGlobal);
    auto targetViewDiff = array::make_view<double, 1>(targetFieldDiffGlobal);

    using for_each = array::helpers::ArrayForEach<0>;

    for_each::apply(std::tie(targetViewDiff ,targetViewCubedSphere, targetViewEqualRegions),
        [&](auto&& diff, auto&& elem1, auto&& elem2) {
        diff = static_cast<double>(elem1) - static_cast<double>(elem2);
        maxDiff = std::max(maxDiff, std::abs(diff));
    });

    targetFunctionSpaceCubedSphere.scatter(targetFieldDiffGlobal, targetFieldDiff);
    targetFieldDiff.haloExchange();

    gmshOutput("interpolation_diff.msh", FieldSet(targetFieldDiff));
    EXPECT_APPROX_EQ(maxDiff, 0., 4 * std::numeric_limits<double>::epsilon());



}



}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
