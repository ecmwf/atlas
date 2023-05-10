/*
 * (C) Copyright 2023 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <numeric>

#include "atlas/array/MakeView.h"
#include "atlas/mesh/BuildMeshFromConnectivities.h"
#include "atlas/mesh/Elements.h"
#include "atlas/mesh/HybridElements.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/Nodes.h"
#include "tests/AtlasTestEnvironment.h"

//#include "atlas/output/Gmsh.h"
//using namespace atlas::output;

namespace atlas {
namespace test {

void test_mesh_setup(const Mesh& mesh, const std::vector<double>& lons, const std::vector<double>& lats,
                     const std::vector<int>& ghosts, const std::vector<gidx_t>& global_indices,
                     const std::vector<idx_t>& remote_indices, const std::vector<int>& partitions,
                     const std::vector<std::array<gidx_t, 3>>& tri_boundary_nodes,
                     const std::vector<gidx_t>& tri_global_indices,
                     const std::vector<std::array<gidx_t, 4>>& quad_boundary_nodes,
                     const std::vector<gidx_t>& quad_global_indices) {
    const auto mesh_xy        = array::make_view<double, 2>(mesh.nodes().xy());
    const auto mesh_lonlat    = array::make_view<double, 2>(mesh.nodes().lonlat());
    const auto mesh_ghost     = array::make_view<int, 1>(mesh.nodes().ghost());
    const auto mesh_gidx      = array::make_view<gidx_t, 1>(mesh.nodes().global_index());
    const auto mesh_ridx      = array::make_view<idx_t, 1>(mesh.nodes().remote_index());
    const auto mesh_partition = array::make_view<int, 1>(mesh.nodes().partition());
    const auto mesh_halo      = array::make_view<int, 1>(mesh.nodes().halo());

    EXPECT(mesh.nodes().size() == lons.size());
    for (size_t i = 0; i < mesh.nodes().size(); ++i) {
        EXPECT(mesh_xy(i, 0) == lons[i]);
        EXPECT(mesh_xy(i, 1) == lats[i]);
        EXPECT(mesh_lonlat(i, 0) == lons[i]);
        EXPECT(mesh_lonlat(i, 1) == lats[i]);
        EXPECT(mesh_ghost(i) == ghosts[i]);
        EXPECT(mesh_gidx(i) == global_indices[i]);
        EXPECT(mesh_ridx(i) == remote_indices[i]);
        EXPECT(mesh_partition(i) == partitions[i]);
        EXPECT(mesh_halo(i) == 0.);
        // Don't expect (or test) any node-to-cell connectivities
    }

    EXPECT(mesh.cells().nb_types() == 2);
    EXPECT(mesh.cells().size() == tri_boundary_nodes.size() + quad_boundary_nodes.size());

    const auto position_of = [&global_indices](const gidx_t idx) {
        const auto& it = std::find(global_indices.begin(), global_indices.end(), idx);
        ATLAS_ASSERT(it != global_indices.end());
        return std::distance(global_indices.begin(), it);
    };

    // Check triangle cell-to-node connectivities
    EXPECT(mesh.cells().elements(0).size() == tri_boundary_nodes.size());
    EXPECT(mesh.cells().elements(0).nb_nodes() == 3);
    for (size_t tri = 0; tri < mesh.cells().elements(0).size(); ++tri) {
        for (size_t node = 0; node < mesh.cells().elements(0).nb_nodes(); ++node) {
            EXPECT(mesh.cells().elements(0).node_connectivity()(tri, node) ==
                   position_of(tri_boundary_nodes[tri][node]));
        }
    }
    // Check quad cell-to-node connectivities
    EXPECT(mesh.cells().elements(1).size() == quad_boundary_nodes.size());
    EXPECT(mesh.cells().elements(1).nb_nodes() == 4);
    for (size_t quad = 0; quad < mesh.cells().elements(1).size(); ++quad) {
        for (size_t node = 0; node < mesh.cells().elements(1).nb_nodes(); ++node) {
            EXPECT(mesh.cells().elements(1).node_connectivity()(quad, node) ==
                   position_of(quad_boundary_nodes[quad][node]));
        }
    }
}

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
    std::vector<idx_t> remote_indices(6);
    std::iota(remote_indices.begin(), remote_indices.end(), 0);  // 0-based numbering
    std::vector<int> partitions(6, 0);                           // all points on proc 0

    // triangles
    std::vector<std::array<gidx_t, 3>> tri_boundary_nodes = {{{3, 6, 5}}, {{3, 4, 6}}};
    std::vector<gidx_t> tri_global_indices                = {1, 2};

    // quads
    std::vector<std::array<gidx_t, 4>> quad_boundary_nodes = {{{1, 2, 3, 5}}};
    std::vector<gidx_t> quad_global_indices                = {3};

    Mesh mesh = build_mesh_from_connectivities(lons, lats, ghosts, global_indices, remote_indices, partitions,
                                               tri_boundary_nodes, tri_global_indices, quad_boundary_nodes,
                                               quad_global_indices);

    test_mesh_setup(mesh, lons, lats, ghosts, global_indices, remote_indices, partitions, tri_boundary_nodes,
                    tri_global_indices, quad_boundary_nodes, quad_global_indices);

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
    std::vector<idx_t> remote_indices(24);
    std::iota(remote_indices.begin(), remote_indices.end(), 0);
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

    Mesh mesh = build_mesh_from_connectivities(lons, lats, ghosts, global_indices, remote_indices, partitions,
                                               tri_boundary_nodes, tri_global_indices, quad_boundary_nodes,
                                               quad_global_indices);

    test_mesh_setup(mesh, lons, lats, ghosts, global_indices, remote_indices, partitions, tri_boundary_nodes,
                    tri_global_indices, quad_boundary_nodes, quad_global_indices);

    //Gmsh gmsh("out.msh", util::Config("coordinates", "xyz"));
    //gmsh.write(mesh);
}

//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
