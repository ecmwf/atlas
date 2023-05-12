/*
 * (C) Copyright 2023 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <algorithm>
#include <array>
#include <vector>

#include "atlas/array/MakeView.h"
#include "atlas/mesh/Elements.h"
#include "atlas/mesh/HybridElements.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/Nodes.h"

namespace atlas {
namespace test {
namespace helper {

void check_mesh_nodes_and_cells(const Mesh& mesh, const std::vector<double>& lons, const std::vector<double>& lats,
                                const std::vector<int>& ghosts, const std::vector<gidx_t>& global_indices,
                                const std::vector<idx_t>& remote_indices, const idx_t remote_index_base,
                                const std::vector<int>& partitions,
                                const std::vector<std::array<gidx_t, 3>>& tri_boundary_nodes,
                                const std::vector<gidx_t>& tri_global_indices,
                                const std::vector<std::array<gidx_t, 4>>& quad_boundary_nodes,
                                const std::vector<gidx_t>& quad_global_indices) {
    const auto mesh_xy        = array::make_view<double, 2>(mesh.nodes().xy());
    const auto mesh_lonlat    = array::make_view<double, 2>(mesh.nodes().lonlat());
    const auto mesh_ghost     = array::make_view<int, 1>(mesh.nodes().ghost());
    const auto mesh_gidx      = array::make_view<gidx_t, 1>(mesh.nodes().global_index());
    const auto mesh_ridx      = array::make_indexview<idx_t, 1>(mesh.nodes().remote_index());
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
        EXPECT(mesh_ridx(i) == remote_indices[i] - remote_index_base);
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

}  // namespace helper
}  // namespace test
}  // namespace atlas
