/*
 * (C) Crown Copyright 2022 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */


#include "atlas/field/FieldSet.h"
#include "atlas/functionspace/CellColumns.h"
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/grid/Grid.h"
#include "atlas/grid/Partitioner.h"
#include "atlas/interpolation/Interpolation.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/meshgenerator/MeshGenerator.h"
#include "atlas/output/Gmsh.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/util/CoordinateEnums.h"
#include "atlas/util/function/VortexRollup.h"

#include "tests/AtlasTestEnvironment.h"


namespace atlas {
namespace test {

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

CASE("cubedsphere_interpolation") {
    // Create a source cubed sphere grid, mesh and functionspace.
    const auto sourceGrid          = Grid("CS-LFR-24");
    const auto sourceMesh          = MeshGenerator("cubedsphere_dual").generate(sourceGrid);
    const auto sourceFunctionspace = functionspace::NodeColumns(sourceMesh);

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

}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
