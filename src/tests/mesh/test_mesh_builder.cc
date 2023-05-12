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
#include "atlas/mesh/Elements.h"
#include "atlas/mesh/HybridElements.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/MeshBuilder.h"
#include "atlas/mesh/Nodes.h"

#include "tests/AtlasTestEnvironment.h"
#include "tests/mesh/helper_mesh_builder.h"

//#include "atlas/output/Gmsh.h"
//using namespace atlas::output;

namespace atlas {
namespace test {

//-----------------------------------------------------------------------------

CASE("test_tiny_mesh") {
    // small regional grid whose cell-centers are connected as (global nodes and cells):
    //
    //   1 - 5 ----- 6
    //   |  3  \ 1 /2|
    //   2 ----- 3 - 4
    //
    std::vector<double> lons{{0.0, 0.0, 10.0, 15.0, 5.0, 15.0}};
    std::vector<double> lats{{5.0, 0.0, 0.0, 0.0, 5.0, 5.0}};

    std::vector<int> ghosts(6, 0);  // all points owned
    std::vector<gidx_t> global_indices(6);
    std::iota(global_indices.begin(), global_indices.end(), 1);  // 1-based numbering
    const idx_t remote_index_base = 0;                           // 0-based numbering
    std::vector<idx_t> remote_indices(6);
    std::iota(remote_indices.begin(), remote_indices.end(), remote_index_base);
    std::vector<int> partitions(6, 0);  // all points on proc 0

    // triangles
    std::vector<std::array<gidx_t, 3>> tri_boundary_nodes = {{{3, 6, 5}}, {{3, 4, 6}}};
    std::vector<gidx_t> tri_global_indices                = {1, 2};

    // quads
    std::vector<std::array<gidx_t, 4>> quad_boundary_nodes = {{{1, 2, 3, 5}}};
    std::vector<gidx_t> quad_global_indices                = {3};

    const MeshBuilder mesh_builder{};
    const Mesh mesh = mesh_builder(lons, lats, ghosts, global_indices, remote_indices, remote_index_base, partitions,
                                   tri_boundary_nodes, tri_global_indices, quad_boundary_nodes, quad_global_indices);

    helper::check_mesh_nodes_and_cells(mesh, lons, lats, ghosts, global_indices, remote_indices, remote_index_base,
                                       partitions, tri_boundary_nodes, tri_global_indices, quad_boundary_nodes,
                                       quad_global_indices);

    //Gmsh gmsh("out.msh", util::Config("coordinates", "xyz"));
    //gmsh.write(mesh);
}

CASE("test_cs_c2_mesh_serial") {
    // coordinates of C2 lfric cubed-sphere grid: grid("CS-LFR-2");
    std::vector<double> lons = {337.5, 22.5,  337.5, 22.5,   // +x
                                67.5,  112.5, 67.5,  112.5,  // +y
                                202.5, 202.5, 157.5, 157.5,  // -x
                                292.5, 292.5, 247.5, 247.5,  // -y
                                315,   45,    225,   135,    // +z
                                315,   225,   45,    135};   // -z

    std::vector<double> lats = {-20.941,  -20.941,  20.941,   20.941,     // +x
                                -20.941,  -20.941,  20.941,   20.941,     // +y
                                -20.941,  20.941,   -20.941,  20.941,     // -x
                                -20.941,  20.941,   -20.941,  20.941,     // -y
                                59.6388,  59.6388,  59.6388,  59.6388,    // +z
                                -59.6388, -59.6388, -59.6388, -59.6388};  // -z

    std::vector<int> ghosts(24, 0);
    std::vector<gidx_t> global_indices(24);
    std::iota(global_indices.begin(), global_indices.end(), 1);
    const idx_t remote_index_base = 1;  // test with 1-based numbering
    std::vector<idx_t> remote_indices(24);
    std::iota(remote_indices.begin(), remote_indices.end(), remote_index_base);
    std::vector<int> partitions(24, 0);

    // triangles
    std::vector<std::array<gidx_t, 3>> tri_boundary_nodes = {//corners
                                                             {{17, 14, 3}},
                                                             {{18, 4, 7}},
                                                             {{20, 8, 12}},
                                                             {{19, 10, 16}},
                                                             {{21, 1, 13}},
                                                             {{23, 5, 2}},
                                                             {{24, 11, 6}},
                                                             {{22, 15, 9}}};
    std::vector<gidx_t> tri_global_indices(8);
    std::iota(tri_global_indices.begin(), tri_global_indices.end(), 1);

    // quads
    std::vector<std::array<gidx_t, 4>> quad_boundary_nodes = {// faces
                                                              {{1, 2, 4, 3}},
                                                              {{5, 6, 8, 7}},
                                                              {{11, 9, 10, 12}},
                                                              {{15, 13, 14, 16}},
                                                              {{17, 18, 20, 19}},
                                                              {{21, 22, 24, 23}},
                                                              // edges between faces
                                                              {{2, 5, 7, 4}},
                                                              {{6, 11, 12, 8}},
                                                              {{9, 15, 16, 10}},
                                                              {{13, 1, 3, 14}},
                                                              {{7, 8, 20, 18}},
                                                              {{12, 10, 19, 20}},
                                                              {{16, 14, 17, 19}},
                                                              {{3, 4, 18, 17}},
                                                              {{23, 24, 6, 5}},
                                                              {{24, 22, 9, 11}},
                                                              {{22, 21, 13, 15}},
                                                              {{21, 23, 2, 1}}};
    std::vector<gidx_t> quad_global_indices(18);
    std::iota(quad_global_indices.begin(), quad_global_indices.end(), 9);  // nb_tris + 1

    const MeshBuilder mesh_builder{};
    const Mesh mesh = mesh_builder(lons, lats, ghosts, global_indices, remote_indices, remote_index_base, partitions,
                                   tri_boundary_nodes, tri_global_indices, quad_boundary_nodes, quad_global_indices);

    helper::check_mesh_nodes_and_cells(mesh, lons, lats, ghosts, global_indices, remote_indices, remote_index_base,
                                       partitions, tri_boundary_nodes, tri_global_indices, quad_boundary_nodes,
                                       quad_global_indices);

    //Gmsh gmsh("out.msh", util::Config("coordinates", "xyz"));
    //gmsh.write(mesh);
}

//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
