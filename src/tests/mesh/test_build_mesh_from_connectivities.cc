/*
 * (C) Copyright 2023 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
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
                     const std::vector<atlas::TriConnectivityData>& tris,
                     const std::vector<atlas::QuadConnectivityData>& quads) {
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
    EXPECT(mesh.cells().size() == tris.size() + quads.size());

    // Check triangle cell-to-node connectivities
    EXPECT(mesh.cells().elements(0).size() == tris.size());
    EXPECT(mesh.cells().elements(0).nb_nodes() == 3);
    for (size_t tri = 0; tri < mesh.cells().elements(0).size(); ++tri) {
        for (size_t node = 0; node < mesh.cells().elements(0).nb_nodes(); ++node) {
            EXPECT(mesh.cells().elements(0).node_connectivity()(tri, node) == tris[tri].boundary_nodes_of_cell[node]);
        }
    }
    // Check quad cell-to-node connectivities
    EXPECT(mesh.cells().elements(1).size() == quads.size());
    EXPECT(mesh.cells().elements(1).nb_nodes() == 4);
    for (size_t quad = 0; quad < mesh.cells().elements(1).size(); ++quad) {
        for (size_t node = 0; node < mesh.cells().elements(1).nb_nodes(); ++node) {
            EXPECT(mesh.cells().elements(1).node_connectivity()(quad, node) ==
                   quads[quad].boundary_nodes_of_cell[node]);
        }
    }
}

//-----------------------------------------------------------------------------

CASE("test_tiny_mesh") {
    // small regional grid whose cell-centers are connected as:
    //
    //   0 - 4 ----- 5
    //   |     \   / |  <-- cells 3,1,2 respectively
    //   1 ----- 2 - 3
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
    std::vector<atlas::TriConnectivityData> tris = {{// local cell index is 0-based; global cell index is 1-based
                                                     {0, 1, {{2, 5, 4}}},
                                                     {1, 2, {{2, 3, 5}}}}};
    // quads
    std::vector<atlas::QuadConnectivityData> quads = {{{2, 3, {{0, 1, 2, 4}}}}};

    Mesh mesh =
        build_mesh_from_connectivities(lons, lats, ghosts, global_indices, remote_indices, partitions, tris, quads);

    test_mesh_setup(mesh, lons, lats, ghosts, global_indices, remote_indices, partitions, tris, quads);

    //Gmsh gmsh("out.msh", util::Config("coordinates", "xyz"));
    //gmsh.write(mesh);
}

CASE("test_cs_c2_mesh_serial") {
    // coordinates of C2 lfric cubed-sphere grid: grid("CS-LFR-2");
    std::vector<double> lons = {{337.5, 22.5,  337.5, 22.5,   // +x
                                 67.5,  112.5, 67.5,  112.5,  // +y
                                 202.5, 202.5, 157.5, 157.5,  // -x
                                 292.5, 292.5, 247.5, 247.5,  // -y
                                 315,   45,    225,   135,    // +z
                                 315,   225,   45,    135}};  // -z

    std::vector<double> lats = {{-20.941,  -20.941,  20.941,   20.941,      // +x
                                 -20.941,  -20.941,  20.941,   20.941,      // +y
                                 -20.941,  20.941,   -20.941,  20.941,      // -x
                                 -20.941,  20.941,   -20.941,  20.941,      // -y
                                 59.6388,  59.6388,  59.6388,  59.6388,     // +z
                                 -59.6388, -59.6388, -59.6388, -59.6388}};  // -z

    std::vector<int> ghosts(24, 0);
    std::vector<gidx_t> global_indices(24);
    std::iota(global_indices.begin(), global_indices.end(), 1);
    std::vector<idx_t> remote_indices(24);
    std::iota(remote_indices.begin(), remote_indices.end(), 0);
    std::vector<int> partitions(24, 0);

    // triangles
    std::vector<atlas::TriConnectivityData> tris = {{// corners
                                                     {0, 1, {{16, 13, 2}}},
                                                     {1, 2, {{17, 3, 6}}},
                                                     {2, 3, {{19, 7, 11}}},
                                                     {3, 4, {{18, 9, 15}}},
                                                     {4, 5, {{20, 0, 12}}},
                                                     {5, 6, {{22, 4, 1}}},
                                                     {6, 7, {{23, 10, 5}}},
                                                     {7, 8, {{21, 14, 8}}}}};
    // quads
    std::vector<atlas::QuadConnectivityData> quads = {{// faces
                                                       {8, 9, {{0, 1, 3, 2}}},
                                                       {9, 10, {{4, 5, 7, 6}}},
                                                       {10, 11, {{10, 8, 9, 11}}},
                                                       {11, 12, {{14, 12, 13, 15}}},
                                                       {12, 13, {{16, 17, 19, 18}}},
                                                       {13, 14, {{20, 21, 23, 22}}},
                                                       // edges between faces
                                                       {14, 15, {{1, 4, 6, 3}}},
                                                       {15, 16, {{5, 10, 11, 7}}},
                                                       {16, 17, {{8, 14, 15, 9}}},
                                                       {17, 18, {{12, 0, 2, 13}}},
                                                       {18, 19, {{6, 7, 19, 17}}},
                                                       {19, 20, {{11, 9, 18, 19}}},
                                                       {20, 21, {{15, 13, 16, 18}}},
                                                       {21, 22, {{2, 3, 17, 16}}},
                                                       {22, 23, {{22, 23, 5, 4}}},
                                                       {23, 24, {{23, 21, 8, 10}}},
                                                       {24, 25, {{21, 20, 12, 14}}},
                                                       {25, 26, {{20, 22, 1, 0}}}}};

    Mesh mesh =
        build_mesh_from_connectivities(lons, lats, ghosts, global_indices, remote_indices, partitions, tris, quads);

    test_mesh_setup(mesh, lons, lats, ghosts, global_indices, remote_indices, partitions, tris, quads);

    //Gmsh gmsh("out.msh", util::Config("coordinates", "xyz"));
    //gmsh.write(mesh);
}

//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
