/*
 * (C) Copyright 2023 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/mesh/BuildMeshFromConnectivities.h"

#include "atlas/array/MakeView.h"
#include "atlas/mesh/ElementType.h"
#include "atlas/mesh/HybridElements.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/util/CoordinateEnums.h"

namespace atlas {

//----------------------------------------------------------------------------------------------------------------------

Mesh build_mesh_from_connectivities(const std::vector<double>& lons, const std::vector<double>& lats,
                                    const std::vector<int>& ghosts, const std::vector<gidx_t>& global_indices,
                                    const std::vector<idx_t>& remote_indices, const std::vector<int>& partitions,
                                    const std::vector<std::array<gidx_t, 3>>& tri_boundary_nodes,
                                    const std::vector<gidx_t>& tri_global_indices,
                                    const std::vector<std::array<gidx_t, 4>>& quad_boundary_nodes,
                                    const std::vector<gidx_t>& quad_global_indices) {
    const size_t nb_points = lons.size();
    ATLAS_ASSERT(nb_points == lats.size());
    ATLAS_ASSERT(nb_points == ghosts.size());
    ATLAS_ASSERT(nb_points == global_indices.size());
    ATLAS_ASSERT(nb_points == remote_indices.size());
    ATLAS_ASSERT(nb_points == partitions.size());

    const size_t nb_tris = tri_boundary_nodes.size();
    ATLAS_ASSERT(nb_tris == tri_global_indices.size());
    const size_t nb_quads = quad_boundary_nodes.size();
    ATLAS_ASSERT(nb_quads == quad_global_indices.size());

    Mesh mesh;

    // nodes

    mesh.nodes().resize(nb_points);
    auto xy        = array::make_view<double, 2>(mesh.nodes().xy());
    auto lonlat    = array::make_view<double, 2>(mesh.nodes().lonlat());
    auto ghost     = array::make_view<int, 1>(mesh.nodes().ghost());
    auto gidx      = array::make_view<gidx_t, 1>(mesh.nodes().global_index());
    auto ridx      = array::make_view<idx_t, 1>(mesh.nodes().remote_index());
    auto partition = array::make_view<int, 1>(mesh.nodes().partition());
    auto halo      = array::make_view<int, 1>(mesh.nodes().halo());

    for (size_t i = 0; i < nb_points; ++i) {
        xy(i, size_t(XX)) = lons[i];
        xy(i, size_t(YY)) = lats[i];
        // identity projection, therefore (lon,lat) = (x,y)
        lonlat(i, size_t(LON)) = lons[i];
        lonlat(i, size_t(LAT)) = lats[i];
        ghost(i)               = ghosts[i];
        gidx(i)                = global_indices[i];
        ridx(i)                = remote_indices[i];
        partition(i)           = partitions[i];
    }
    halo.assign(0);

    // cells

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
    const auto position_of = [&global_indices](const gidx_t idx) {
        const auto& it = std::find(global_indices.begin(), global_indices.end(), idx);
        ATLAS_ASSERT(it != global_indices.end());
        return std::distance(global_indices.begin(), it);
    };
    // Find array of positions inside global_indices from an array of ixds
    // (In C++20 we can template the lambda on array size; for now we use an auto argument...)
    const auto positions_of = [&position_of](const auto idxs) {
        std::array<idx_t, idxs.size()> result{};  // go from input gidx_t to output idx_t
        std::transform(idxs.begin(), idxs.end(), result.begin(), position_of);
        return result;
    };

    size_t idx = 0;
    if (add_tris) {
        for (size_t tri = 0; tri < nb_tris; ++tri) {
            node_connectivity.set(idx, positions_of(tri_boundary_nodes[tri]).data());
            cells_gidx(idx) = tri_global_indices[tri];
            idx++;
        }
    }
    if (add_quads) {
        for (size_t quad = 0; quad < nb_quads; ++quad) {
            node_connectivity.set(idx, positions_of(quad_boundary_nodes[quad]).data());
            cells_gidx(idx) = quad_global_indices[quad];
            idx++;
        }
    }

    ATLAS_ASSERT(idx == nb_tris + nb_quads);

    cells_part.assign(atlas::mpi::comm().rank());

    return mesh;
}

//----------------------------------------------------------------------------------------------------------------------

}  // namespace atlas
