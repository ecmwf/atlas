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
#include "atlas/functionspace/PointCloud.h"
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
#include "atlas/util/CoordinateEnums.h"
#include "atlas/util/function/VortexRollup.h"

#include "tests/AtlasTestEnvironment.h"


namespace atlas {
namespace test {

struct CubedSphereInterpolationFixture {
    atlas::Grid sourceGrid_ = Grid("CS-LFR-24");
    atlas::Mesh sourceMesh_ = MeshGenerator("cubedsphere_dual").generate(sourceGrid_);
    atlas::FunctionSpace sourceFunctionSpace_ = functionspace::NodeColumns(sourceMesh_);
    atlas::grid::Partitioner targetPartitioner_ =
        grid::MatchingPartitioner(sourceMesh_, util::Config("type", "cubedsphere"));
    atlas::Grid targetGrid_ = Grid("O24");
    atlas::Mesh targetMesh_ = MeshGenerator("structured").generate(targetGrid_, targetPartitioner_);
    atlas::FunctionSpace targetFunctionSpace_ = functionspace::NodeColumns(targetMesh_);
};

void gmshOutput(const std::string& fileName, const FieldSet& fieldSet) {


    const auto functionSpace = fieldSet[0].functionspace();
    const auto mesh = functionspace::NodeColumns(functionSpace).mesh();

    const auto gmshConfig =
        util::Config("coordinates", "xyz") | util::Config("ghost", true) | util::Config("info", true);
    const auto gmsh = output::Gmsh(fileName, gmshConfig);
    gmsh.write(mesh);
    gmsh.write(fieldSet, functionSpace);
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

CASE("temporary const lvalue referene test") {

    const auto grid = CubedSphereInterpolationFixture{}.sourceGrid_;

    auto iterator = grid.lonlat().begin();
    const auto begin = iterator;
    for (idx_t i = 0; i < grid.size(); ++i) {

        // This does not dangle, because C++ makes no sense!
        const auto& lonLat = *(begin + i);
        const auto refLonLat = *iterator;
        ++iterator;

        EXPECT_APPROX_EQ(lonLat.lon(), refLonLat.lon(), 1e-14);
        EXPECT_APPROX_EQ(lonLat.lat(), refLonLat.lat(), 1e-14);

        // Only now, *(begin + i), do you have my permission to die.
    }


}

CASE("cubedsphere_to_cubedsphere_interpolation") {

    // Check that interpolation scheme correctly maps a cubedsphere on to itself.

    const auto fixture = CubedSphereInterpolationFixture{};
    auto sourceField = fixture.sourceFunctionSpace_.createField<double>(option::name("test_field"));
    {
        const auto lonlat = array::make_view<double, 2>(fixture.sourceFunctionSpace_.lonlat());
        auto view         = array::make_view<double, 1>(sourceField);
        for (idx_t i = 0; i < fixture.sourceFunctionSpace_.size(); ++i) {
            view(i) = util::function::vortex_rollup(lonlat(i, LON), lonlat(i, LAT), 1.);
        }
    }
    const auto targetMesh = MeshGenerator("cubedsphere_dual").generate(fixture.sourceGrid_, fixture.targetPartitioner_);
    const auto targetFunctionSpace = functionspace::CubedSphereNodeColumns(targetMesh);

    const auto scheme = util::Config("type", "cubedsphere-bilinear");
    const auto interp = Interpolation(scheme, fixture.sourceFunctionSpace_, targetFunctionSpace);

    auto targetField = targetFunctionSpace.createField<double>(option::name("test_field"));
    interp.execute(sourceField, targetField);
    targetField.haloExchange();

    // iterate over target field and check that error is zero.
    {
        const auto lonlat = array::make_view<double, 2>(targetFunctionSpace.lonlat());
        auto targetView   = array::make_view<double, 1>(targetField);
        for (idx_t i = 0; i < targetFunctionSpace.size(); ++i) {
            const auto val = util::function::vortex_rollup(lonlat(i, LON), lonlat(i, LAT), 1.);
            EXPECT_APPROX_EQ(targetView(i), val, 1e-14);
        }
    }

    // output source and target fields.
    gmshOutput("cubedsphere_to_cubedsphere_source.msh", FieldSet{sourceField});
    gmshOutput("cubedsphere_to_cubedsphere_target.msh", FieldSet{targetField});
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
    const auto scheme = util::Config("type", "cubedsphere-bilinear") | util::Config("adjoint", true);
    const auto interp = Interpolation(scheme, fixture.sourceFunctionSpace_, fixture.targetFunctionSpace_);

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

}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
