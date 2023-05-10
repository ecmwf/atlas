/*
 * (C) Copyright 2023 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <algorithm>
#include <array>
#include <iterator>
#include <type_traits>
#include <vector>

#include "atlas/array/MakeView.h"
#include "atlas/mesh/BuildMeshFromConnectivities.h"
#include "atlas/mesh/Elements.h"
#include "atlas/mesh/HybridElements.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/parallel/mpi/mpi.h"
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

CASE("test_cs_c2_mesh_parallel") {
    ATLAS_ASSERT(mpi::comm().size() == 6);
    const int rank = mpi::comm().rank();

    // Coordinates of the C2 LFRic cubed-sphere grid: grid("CS-LFR-2");
    const std::vector<double> global_lons = {337.5, 22.5,  337.5, 22.5,   // +x
                                             67.5,  112.5, 67.5,  112.5,  // +y
                                             202.5, 202.5, 157.5, 157.5,  // -x
                                             292.5, 292.5, 247.5, 247.5,  // -y
                                             315,   45,    225,   135,    // +z
                                             315,   225,   45,    135};   // -z

    const std::vector<double> global_lats = {-20.941,  -20.941,  20.941,   20.941,     // +x
                                             -20.941,  -20.941,  20.941,   20.941,     // +y
                                             -20.941,  20.941,   -20.941,  20.941,     // -x
                                             -20.941,  20.941,   -20.941,  20.941,     // -y
                                             59.6388,  59.6388,  59.6388,  59.6388,    // +z
                                             -59.6388, -59.6388, -59.6388, -59.6388};  // -z

    const std::vector<int> global_partitions = {0, 0, 0, 0,   // +x
                                                1, 1, 1, 1,   // +y
                                                2, 2, 2, 2,   // -x
                                                3, 3, 3, 3,   // -y
                                                4, 4, 4, 4,   // +z
                                                5, 5, 5, 5};  // -z

    const std::vector<int> local_indices = {0, 1, 2, 3,   // +x
                                            0, 1, 2, 3,   // +y
                                            0, 1, 2, 3,   // -x
                                            0, 1, 2, 3,   // -y
                                            0, 1, 2, 3,   // +z
                                            0, 1, 2, 3};  // -z

    std::vector<gidx_t> global_indices;
    std::vector<std::array<gidx_t, 3>> tri_boundary_nodes{};
    std::vector<gidx_t> tri_global_indices{};
    std::vector<std::array<gidx_t, 4>> quad_boundary_nodes{};
    std::vector<gidx_t> quad_global_indices{};

    if (rank == 0) {
        // global indices for points and cells are 1-based
        global_indices = {1, 2, 3, 4, 5, 7, 17, 18};
        tri_boundary_nodes.push_back({{18, 4, 7}});
        tri_global_indices = {2};
        quad_boundary_nodes.push_back({{1, 2, 4, 3}});
        quad_boundary_nodes.push_back({{2, 5, 7, 4}});
        quad_boundary_nodes.push_back({{3, 4, 18, 17}});
        quad_global_indices = {9, 15, 22};
    }
    else if (rank == 1) {
        global_indices = {5, 6, 7, 8, 11, 12, 18, 20};
        tri_boundary_nodes.push_back({{20, 8, 12}});
        tri_global_indices = {3};
        quad_boundary_nodes.push_back({{5, 6, 8, 7}});
        quad_boundary_nodes.push_back({{6, 11, 12, 8}});
        quad_boundary_nodes.push_back({{7, 8, 20, 18}});
        quad_global_indices = {10, 16, 19};
    }
    else if (rank == 2) {
        global_indices = {9, 10, 11, 12, 15, 16, 22, 24};
        tri_boundary_nodes.push_back({{22, 15, 9}});
        tri_global_indices = {8};
        quad_boundary_nodes.push_back({{11, 9, 10, 12}});
        quad_boundary_nodes.push_back({{9, 15, 16, 10}});
        quad_boundary_nodes.push_back({{24, 22, 9, 11}});
        quad_global_indices = {11, 17, 24};
    }
    else if (rank == 3) {
        global_indices = {1, 3, 13, 14, 15, 16, 21, 22};
        tri_boundary_nodes.push_back({{21, 1, 13}});
        tri_global_indices = {5};
        quad_boundary_nodes.push_back({{15, 13, 14, 16}});
        quad_boundary_nodes.push_back({{13, 1, 3, 14}});
        quad_boundary_nodes.push_back({{22, 21, 13, 15}});
        quad_global_indices = {12, 18, 25};
    }
    else if (rank == 4) {
        global_indices = {3, 10, 12, 14, 16, 17, 18, 19, 20};
        tri_boundary_nodes.push_back({{17, 14, 3}});
        tri_boundary_nodes.push_back({{19, 10, 16}});
        tri_global_indices = {1, 4};
        quad_boundary_nodes.push_back({{17, 18, 20, 19}});
        quad_boundary_nodes.push_back({{12, 10, 19, 20}});
        quad_boundary_nodes.push_back({{16, 14, 17, 19}});
        quad_global_indices = {13, 20, 21};
    }
    else {  // rank == 5
        global_indices = {1, 2, 5, 6, 11, 21, 22, 23, 24};
        tri_boundary_nodes.push_back({{23, 5, 2}});
        tri_boundary_nodes.push_back({{24, 11, 6}});
        tri_global_indices = {6, 7};
        quad_boundary_nodes.push_back({{21, 22, 24, 23}});
        quad_boundary_nodes.push_back({{23, 24, 6, 5}});
        quad_boundary_nodes.push_back({{21, 23, 2, 1}});
        quad_global_indices = {14, 23, 26};
    }

    // Compute (local subset of) {lons,lats,ghosts,partitions} from (local subset of) global indices
    std::vector<double> lons;
    std::transform(global_indices.begin(), global_indices.end(), std::back_inserter(lons),
                   [&global_lons](const gidx_t idx) { return global_lons[idx - 1]; });

    std::vector<double> lats;
    std::transform(global_indices.begin(), global_indices.end(), std::back_inserter(lats),
                   [&global_lats](const gidx_t idx) { return global_lats[idx - 1]; });

    std::vector<int> ghosts;
    std::transform(
        global_indices.begin(), global_indices.end(), std::back_inserter(ghosts),
        [&global_partitions, &rank](const gidx_t idx) { return static_cast<int>(global_partitions[idx - 1] != rank); });

    std::vector<idx_t> remote_indices;
    std::transform(global_indices.begin(), global_indices.end(), std::back_inserter(remote_indices),
                   [&local_indices](const gidx_t idx) { return static_cast<idx_t>(local_indices[idx - 1]); });

    std::vector<int> partitions;
    std::transform(global_indices.begin(), global_indices.end(), std::back_inserter(partitions),
                   [&global_partitions](const gidx_t idx) { return global_partitions[idx - 1]; });

    Mesh mesh = build_mesh_from_connectivities(lons, lats, ghosts, global_indices, remote_indices, partitions,
                                               tri_boundary_nodes, tri_global_indices, quad_boundary_nodes,
                                               quad_global_indices);

    test_mesh_setup(mesh, lons, lats, ghosts, global_indices, remote_indices, partitions, tri_boundary_nodes,
                    tri_global_indices, quad_boundary_nodes, quad_global_indices);

    //Gmsh gmsh("fh_debug_out.msh", util::Config("coordinates", "xyz"));
    //gmsh.write(mesh);
}

//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
