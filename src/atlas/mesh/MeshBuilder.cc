/*
 * (C) Copyright 2023 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "atlas/mesh/MeshBuilder.h"

#include <algorithm>

#include "atlas/array/MakeView.h"
#include "atlas/grid/Grid.h"
#include "atlas/grid/Iterator.h"
#include "atlas/grid/StructuredGrid.h"
#include "atlas/grid/UnstructuredGrid.h"
#include "atlas/mesh/ElementType.h"
#include "atlas/mesh/HybridElements.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/util/CoordinateEnums.h"

namespace atlas {

namespace detail {
atlas::UnstructuredGrid assemble_unstructured_grid(size_t nb_nodes, const double lons[], const double lats[],
                                                   const int ghosts[], const eckit::mpi::Comm& comm) {

    // First serialize owned lons and lats into single vector
    const size_t nb_owned_nodes = std::count(ghosts, ghosts + nb_nodes, 0);
    std::vector<double> owned_lonlats(2 * nb_owned_nodes);
    int counter = 0;
    for (size_t i = 0; i < nb_nodes; ++i) {
        if (ghosts[i] == 0) {
            owned_lonlats[counter] = lons[i];
            counter++;
            owned_lonlats[counter] = lats[i];
            counter++;
        }
    }
    ATLAS_ASSERT(counter == 2 * nb_owned_nodes);

    // Gather points across MPI ranks
    size_t nb_nodes_global = 0;
    comm.allReduce(nb_owned_nodes, nb_nodes_global, eckit::mpi::sum());

    std::vector<double> global_lonlats(2 * nb_nodes_global);
    eckit::mpi::Buffer<double> buffer(comm.size());
    comm.allGatherv(owned_lonlats.begin(), owned_lonlats.end(), buffer);
    global_lonlats = std::move(buffer.buffer);

    std::vector<atlas::PointXY> points(nb_nodes_global);
    for (size_t i = 0; i < nb_nodes_global; ++i) {
        points[i] = atlas::PointXY({global_lonlats[2*i], global_lonlats[2*i + 1]});
    }

    return atlas::UnstructuredGrid(new std::vector<atlas::PointXY>(points.begin(), points.end()));
}

// Check if grid points match mesh points -- not obviously true as the MeshBuilder sets the Grid
// independently from the Mesh. This check can give confidence that the grid and mesh are specified
// consistently with each other.
//
// This check makes a few assumptions
// - the global_indices passed to the MeshBuilder are a 1-based contiguous index
// - the Grid is one of StructuredGrid or UnstructedGrid. This restriction is because a
//   CubedSphereGrid doesn't present a simple interface for the coordinates at the N'th point.
// - the Mesh grid points are the grid nodes, and not the grid cell-centers (e.g., HEALPix or CubedSphere meshes)
void validate_grid_vs_mesh(const atlas::Grid& grid, size_t nb_nodes, const double lons[], const double lats[],
                           const int ghosts[], const gidx_t global_indices[], const eckit::mpi::Comm& comm) {
    // Check assumption that global_indices look like a 1-based contiguous index
    const size_t nb_owned_nodes = std::count(ghosts, ghosts + nb_nodes, 0);
    size_t nb_nodes_global = 0;
    comm.allReduce(nb_owned_nodes, nb_nodes_global, eckit::mpi::sum());
    for (size_t i = 0; i < nb_nodes; ++i) {
        if (ghosts[i] == 0) {
            // Check global_indices is consistent with a 1-based contiguous index over nodes
            ATLAS_ASSERT(global_indices[i] >= 1);
            ATLAS_ASSERT(global_indices[i] <= nb_nodes_global);
        }
    }

    double lonlat[2];
    const auto equal_within_roundoff = [](const double a, const double b) -> bool {
        return fabs(a - b) <= 360.0 * 1.0e-16;
    };

    // Check lonlats for each supported grid type
    if (grid.type() == "unstructured") {
        const atlas::UnstructuredGrid ugrid(grid);
        for (size_t i = 0; i < nb_nodes; ++i) {
            if (ghosts[i] == 0) {
                ugrid.lonlat(global_indices[i] - 1, lonlat);
                if (!equal_within_roundoff(lonlat[0], lons[i]) || !equal_within_roundoff(lonlat[1], lats[i])) {
                    throw_Exception("In MeshBuilder: UnstructuredGrid from config does not match mesh coordinates", Here());
                }
            }
        }
    } else if (grid.type() == "structured") {
        const atlas::StructuredGrid sgrid(grid);
        for (size_t i = 0; i < nb_nodes; ++i) {
            if (ghosts[i] == 0) {
                idx_t i, j;
                sgrid.index2ij(global_indices[i] - 1, i, j);
                sgrid.lonlat(i, j, lonlat);
                if (!equal_within_roundoff(lonlat[0], lons[i]) || !equal_within_roundoff(lonlat[1], lats[i])) {
                    throw_Exception("In MeshBuilder: StructuredGrid from config does not match mesh coordinates", Here());
                }
            }
        }
    } else {
        for (size_t i = 0; i < nb_nodes; ++i) {
            if (ghosts[i] == 0) {
                auto point = *(grid.lonlat().begin() + static_cast<size_t>(global_indices[i] - 1));
                if (!equal_within_roundoff(point.lon(), lons[i]) || !equal_within_roundoff(point.lat(), lats[i])) {
                    throw_Exception("In MeshBuilder: Grid from config does not match mesh coordinates", Here());
                }
            }
        }
    }
}
}  // namespace detail

namespace mesh {

//----------------------------------------------------------------------------------------------------------------------

Mesh MeshBuilder::operator()(const std::vector<double>& lons, const std::vector<double>& lats,
                             const std::vector<int>& ghosts, const std::vector<gidx_t>& global_indices,
                             const std::vector<idx_t>& remote_indices, const idx_t remote_index_base,
                             const std::vector<int>& partitions,
                             const std::vector<std::array<gidx_t, 3>>& tri_boundary_nodes,
                             const std::vector<gidx_t>& tri_global_indices,
                             const std::vector<std::array<gidx_t, 4>>& quad_boundary_nodes,
                             const std::vector<gidx_t>& quad_global_indices,
                             const eckit::Configuration& config) const {
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
                      reinterpret_cast<const gidx_t*>(quad_boundary_nodes.data()), quad_global_indices.data(),
                      config);
}

Mesh MeshBuilder::operator()(size_t nb_nodes, const double lons[], const double lats[], const int ghosts[],
                             const gidx_t global_indices[], const idx_t remote_indices[], const idx_t remote_index_base,
                             const int partitions[], size_t nb_tris, const gidx_t tri_boundary_nodes[],
                             const gidx_t tri_global_indices[], size_t nb_quads, const gidx_t quad_boundary_nodes[],
                             const gidx_t quad_global_indices[],
                             const eckit::Configuration& config) const {
    // Get MPI comm from config name or fall back to atlas default comm
    auto mpi_comm_name = [](const auto& config) {
        return config.getString("mpi_comm", atlas::mpi::comm().name());
    };
    const eckit::mpi::Comm& comm = eckit::mpi::comm(mpi_comm_name(config).c_str());

    Mesh mesh{};

    // Setup a grid, if requested via config argument
    if (config.has("grid")) {
        atlas::Grid grid;
        if (config.has("grid.type") && (config.getString("grid.type") == "unstructured")
            && !config.has("grid.xy")) {
            // Assemble the unstructured grid by gathering input lons,lats across ranks
            grid = ::atlas::detail::assemble_unstructured_grid(nb_nodes, lons, lats, ghosts, comm);
        } else {
            // Build grid directly from config
            grid = atlas::Grid(config.getSubConfiguration("grid"));
            const bool validate = config.getBool("validate", false);
            if (validate) {
                ::atlas::detail::validate_grid_vs_mesh(grid, nb_nodes, lons, lats, ghosts, global_indices, comm);
            }
        }
        mesh.setGrid(grid);
    }

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
    comm.allReduce(nb_tris, sum_nb_tris, eckit::mpi::sum());
    const bool add_tris = (sum_nb_tris > 0);

    size_t sum_nb_quads = 0;
    comm.allReduce(nb_quads, sum_nb_quads, eckit::mpi::sum());
    const bool add_quads = (sum_nb_quads > 0);

    if (add_tris) {
        mesh.cells().add(mesh::ElementType::create("Triangle"), nb_tris);
    }
    if (add_quads) {
        mesh.cells().add(mesh::ElementType::create("Quadrilateral"), nb_quads);
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

    cells_part.assign(comm.rank());

    return mesh;
}

//----------------------------------------------------------------------------------------------------------------------

}  // namespace mesh
}  // namespace atlas
