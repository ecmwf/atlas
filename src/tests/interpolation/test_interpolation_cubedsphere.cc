/*
 * (C) Crown Copyright 2022 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */


#include "atlas/field/FieldSet.h"
#include "atlas/functionspace/CellColumns.h"
#include "atlas/grid/Grid.h"
#include "atlas/grid/Partitioner.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/meshgenerator/MeshGenerator.h"
#include "atlas/output/Gmsh.h"
#include "atlas/redistribution/Redistribution.h"

#include "tests/AtlasTestEnvironment.h"





namespace atlas {
namespace test {

CASE("cubedsphere_matching_mesh") {

    // Set grid, mesh and functionspace.
    const auto sourceGrid = Grid("CS-LFR-C-24");
    const auto sourceMesh = MeshGenerator("cubedsphere_dual").generate(sourceGrid);

    const auto partitioner = grid::MatchingPartitioner(sourceMesh, util::Config("type", "cubedsphere"));
    const auto targetGrid = Grid("O32");


    const auto interMesh = MeshGenerator("structured").generate(targetGrid, partitioner);
    const auto interFunctionspace = functionspace::NodeColumns(interMesh);

    auto interField = interFunctionspace.createField<int>(option::name("partition"));
    {
        auto view = array::make_view<int, 1>(interField);
        for (idx_t i = 0; i < interFunctionspace.size(); ++i) {
            view(i) = mpi::rank();
        }
    }
    interField.haloExchange();

    const auto targetMesh = MeshGenerator("structured").generate(targetGrid);
    const auto targetFunctionspace = functionspace::NodeColumns(targetMesh);
    auto targetField = targetFunctionspace.createField<int>(option::name("partition"));

    const auto redist = Redistribution(interFunctionspace, targetFunctionspace);
    redist.execute(interField, targetField);
    targetField.haloExchange();

    const auto gmshConfig =
        util::Config("coordinates", "xyz") | util::Config("ghost", true) | util::Config("info", true);
    const auto intermediateGmsh = output::Gmsh("intermediate.msh", gmshConfig);
    intermediateGmsh.write(interMesh);
    intermediateGmsh.write(interField, interFunctionspace);

    const auto targetGmsh = output::Gmsh("target.msh", gmshConfig);
    targetGmsh.write(targetMesh);
    targetGmsh.write(targetField, targetFunctionspace);



}

}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
