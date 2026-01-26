/*
 * (C) Copyright 2023 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <array>
#include <numeric>
#include <vector>

#include "atlas/array/MakeView.h"
#include "atlas/grid.h"
#include "atlas/mesh/Elements.h"
#include "atlas/mesh/HybridElements.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/MeshBuilder.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/util/Config.h"
#include "atlas/output/Gmsh.h"

#include "tests/AtlasTestEnvironment.h"

using namespace atlas::mesh;

//#include "atlas/output/Gmsh.h"
//using namespace atlas::output;

namespace atlas {
namespace test {

//-----------------------------------------------------------------------------

CASE("test_tiny_mesh") {
    // small regional grid whose cell-centers are connected as (global nodes and cells):
    //
    //   1 ---- 5 ----- 6
    //   | 3 / 4 | 1 /2 |
    //   2 ----- 3 ---- 4
    //
    gidx_t global_index_base = 1;
    size_t nb_nodes = 6;
    std::vector<double> lon{{0.0, 0.0, 10.0, 15.0, 5.0, 15.0}};
    std::vector<double> lat{{5.0, 0.0, 0.0, 0.0, 5.0, 5.0}};
    std::vector<double> x(nb_nodes);
    std::vector<double> y(nb_nodes);
    for (size_t j=0; j<nb_nodes; ++j) {
        x[j] = lon[j] / 10.;
        y[j] = lat[j] / 10.;
    }
    std::vector<gidx_t> global_index(6);
    std::iota(global_index.begin(), global_index.end(), global_index_base);

    // triangles
    size_t nb_triags = 4;
    std::vector<std::array<gidx_t, 3>> triag_nodes_global = {{{3, 6, 5}}, {{3, 4, 6}}, {{2, 5, 1}}, {{2, 3, 5}}};
    std::vector<gidx_t> triag_global_index                = {1, 2, 3, 4};

    const TriangularMeshBuilder mesh_builder{};
    constexpr size_t stride_1 = 1;
    const Mesh mesh = mesh_builder(nb_nodes, global_index.data(),
                                   x.data(), y.data(), stride_1, stride_1, lon.data(), lat.data(), stride_1, stride_1,
                                   nb_triags, triag_global_index.data(), triag_nodes_global.data()->data(),
                                   global_index_base);

    output::Gmsh gmsh("out.msh", util::Config("coordinates", "xy"));
    gmsh.write(mesh);
}

//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
