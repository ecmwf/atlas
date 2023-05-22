/*
 * (C) Copyright 2023 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "atlas/mesh/MeshBuilder.h"

#include "atlas/array/MakeView.h"
#include "atlas/mesh/ElementType.h"
#include "atlas/mesh/HybridElements.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/util/CoordinateEnums.h"

namespace atlas {
namespace mesh {

//----------------------------------------------------------------------------------------------------------------------

Mesh MeshBuilder::operator()(const std::vector<double>& lons, const std::vector<double>& lats,
                             const std::vector<int>& ghosts, const std::vector<gidx_t>& global_indices,
                             const std::vector<idx_t>& remote_indices, const idx_t remote_index_base,
                             const std::vector<int>& partitions,
                             const std::vector<std::array<gidx_t, 3>>& tri_boundary_nodes,
                             const std::vector<gidx_t>& tri_global_indices,
                             const std::vector<std::array<gidx_t, 4>>& quad_boundary_nodes,
                             const std::vector<gidx_t>& quad_global_indices) const {
    const size_t nb_nodes = global_indices.size();
    const size_t nb_tris  = tri_global_indices.size();
    const size_t nb_quads = quad_global_indices.size();

    ATLAS_ASSERT(nb_nodes == lons.size());
    ATLAS_ASSERT(nb_nodes == lats.size());
    ATLAS_ASSERT(nb_nodes == ghosts.size());
    ATLAS_ASSERT(nb_nodes == remote_indices.size());
    ATLAS_ASSERT(nb_nodes == partitions.size());
    ATLAS_ASSERT(nb_tris == tri_boundary_nodes.size());
    ATLAS_ASSERT(nb_quads == quad_boundary_nodes.size());

    return operator()(nb_nodes, lons.data(), lats.data(), ghosts.data(), global_indices.data(), remote_indices.data(),
                      remote_index_base, partitions.data(), nb_tris,
                      reinterpret_cast<const gidx_t*>(tri_boundary_nodes.data()), tri_global_indices.data(), nb_quads,
                      reinterpret_cast<const gidx_t*>(quad_boundary_nodes.data()), quad_global_indices.data());
}

Mesh MeshBuilder::operator()(size_t nb_nodes, const double lons[], const double lats[], const int ghosts[],
                             const gidx_t global_indices[], const idx_t remote_indices[], const idx_t remote_index_base,
                             const int partitions[], size_t nb_tris, const gidx_t tri_boundary_nodes[],
                             const gidx_t tri_global_indices[], size_t nb_quads, const gidx_t quad_boundary_nodes[],
                             const gidx_t quad_global_indices[]) const {
    Mesh mesh{};

    // Populate node data

    mesh.nodes().resize(nb_nodes);
    auto xy        = array::make_view<double, 2>(mesh.nodes().xy());
    auto lonlat    = array::make_view<double, 2>(mesh.nodes().lonlat());
    auto ghost     = array::make_view<int, 1>(mesh.nodes().ghost());
    auto gidx      = array::make_view<gidx_t, 1>(mesh.nodes().global_index());
    auto ridx      = array::make_indexview<idx_t, 1>(mesh.nodes().remote_index());
    auto partition = array::make_view<int, 1>(mesh.nodes().partition());
    auto halo      = array::make_view<int, 1>(mesh.nodes().halo());

    for (size_t i = 0; i < nb_nodes; ++i) {
        xy(i, size_t(XX)) = lons[i];
        xy(i, size_t(YY)) = lats[i];
        // Identity projection, therefore (lon,lat) = (x,y)
        lonlat(i, size_t(LON)) = lons[i];
        lonlat(i, size_t(LAT)) = lats[i];
        ghost(i)               = ghosts[i];
        gidx(i)                = global_indices[i];
        ridx(i)                = remote_indices[i] - remote_index_base;
        partition(i)           = partitions[i];
    }
    halo.assign(0);

    // Populate cell/element data

    // First, count how many cells of each type are on this processor
    // Then optimize away the element type if globally nb_tris or nb_quads is zero
    size_t sum_nb_tris = 0;
    atlas::mpi::comm().allReduce(nb_tris, sum_nb_tris, eckit::mpi::sum());
    const bool add_tris = (sum_nb_tris > 0);

    size_t sum_nb_quads = 0;
    atlas::mpi::comm().allReduce(nb_quads, sum_nb_quads, eckit::mpi::sum());
    const bool add_quads = (sum_nb_quads > 0);

    if (add_tris) {
        mesh.cells().add(new mesh::temporary::Triangle(), nb_tris);
    }
    if (add_quads) {
        mesh.cells().add(new mesh::temporary::Quadrilateral(), nb_quads);
    }

    atlas::mesh::HybridElements::Connectivity& node_connectivity = mesh.cells().node_connectivity();
    auto cells_part = array::make_view<int, 1>(mesh.cells().partition());
    auto cells_gidx = array::make_view<gidx_t, 1>(mesh.cells().global_index());

    // Find position of idx inside global_indices
    const auto position_of = [&nb_nodes, &global_indices](const gidx_t idx) {
        const auto& it = std::find(global_indices, global_indices + nb_nodes, idx);
        ATLAS_ASSERT(it != global_indices + nb_nodes);
        return std::distance(global_indices, it);
    };

    size_t idx = 0;
    if (add_tris) {
        idx_t buffer[3];
        for (size_t tri = 0; tri < nb_tris; ++tri) {
            for (size_t i = 0; i < 3; ++i) {
                buffer[i] = position_of(tri_boundary_nodes[3 * tri + i]);
            }
            node_connectivity.set(idx, buffer);
            cells_gidx(idx) = tri_global_indices[tri];
            idx++;
        }
    }
    if (add_quads) {
        idx_t buffer[4];
        for (size_t quad = 0; quad < nb_quads; ++quad) {
            for (size_t i = 0; i < 4; ++i) {
                buffer[i] = position_of(quad_boundary_nodes[4 * quad + i]);
            }
            node_connectivity.set(idx, buffer);
            cells_gidx(idx) = quad_global_indices[quad];
            idx++;
        }
    }

    ATLAS_ASSERT(idx == nb_tris + nb_quads);

    cells_part.assign(atlas::mpi::comm().rank());

    return mesh;
}

//----------------------------------------------------------------------------------------------------------------------

}  // namespace mesh
}  // namespace atlas
