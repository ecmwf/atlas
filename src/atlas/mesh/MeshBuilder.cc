/*
 * (C) Copyright 2023 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "atlas/mesh/MeshBuilder.h"

#include <algorithm>
#include <numeric>

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
atlas::UnstructuredGrid assemble_unstructured_grid(size_t nb_nodes, const double lon[], const double lat[],
                                                   const int ghost[], const eckit::mpi::Comm& comm) {

    // First serialize owned lons and lats into single vector
    const size_t nb_owned_nodes = std::count(ghost, ghost + nb_nodes, 0);
    std::vector<double> owned_lonlat(2 * nb_owned_nodes);
    int counter = 0;
    for (size_t n = 0; n < nb_nodes; ++n) {
        if (ghost[n] == 0) {
            owned_lonlat[counter] = lon[n];
            counter++;
            owned_lonlat[counter] = lat[n];
            counter++;
        }
    }
    ATLAS_ASSERT(counter == 2 * nb_owned_nodes);

    // Gather points across MPI ranks
    size_t nb_nodes_global = 0;
    comm.allReduce(nb_owned_nodes, nb_nodes_global, eckit::mpi::sum());

    std::vector<double> global_lonlat(2 * nb_nodes_global);
    eckit::mpi::Buffer<double> buffer(comm.size());
    comm.allGatherv(owned_lonlat.begin(), owned_lonlat.end(), buffer);
    global_lonlat = std::move(buffer.buffer);

    std::vector<atlas::PointXY> points(nb_nodes_global);
    for (size_t n = 0; n < nb_nodes_global; ++n) {
        points[n] = atlas::PointXY({global_lonlat[2*n], global_lonlat[2*n + 1]});
    }

    return atlas::UnstructuredGrid(new std::vector<atlas::PointXY>(points.begin(), points.end()));
}

// Check if grid points match mesh points -- not obviously true as the MeshBuilder sets the Grid
// independently from the Mesh. This check can give confidence that the grid and mesh are specified
// consistently with each other.
//
// This check makes a few assumptions
// - the global_index passed to the MeshBuilder are a global_index_base-based contiguous index
// - the Grid is one of StructuredGrid or UnstructedGrid. This restriction is because a
//   CubedSphereGrid doesn't present a simple interface for the coordinates at the N'th point.
// - the Mesh grid points are the grid nodes, and not the grid cell-centers (e.g., HEALPix or CubedSphere meshes)
void validate_grid_vs_mesh(const eckit::mpi::Comm& comm, const atlas::Grid& grid,
                           size_t nb_nodes, const double lon[], const double lat[], size_t lonstride, size_t latstride,
                           const int ghost[], const gidx_t global_index[], gidx_t global_index_base) {
    // Check assumption that global_index look like a global_index_based contiguous index
    const size_t nb_owned_nodes = std::count(ghost, ghost + nb_nodes, 0);
    size_t nb_nodes_global = 0;
    comm.allReduce(nb_owned_nodes, nb_nodes_global, eckit::mpi::sum());
    for (size_t n = 0; n < nb_nodes; ++n) {
        if (ghost[n] == 0) {
            // Check global_index is consistent with a global_index_base contiguous index over nodes
            auto g = global_index[n] - global_index_base;
            ATLAS_ASSERT(g >= 0);
            ATLAS_ASSERT(g < nb_nodes_global);
        }
    }

    double lonlat[2];
    const auto equal_within_roundoff = [](const double a, const double b) -> bool {
        return std::abs(a - b) <= 360.0 * 1.0e-16;
    };

    // Check lonlats for each supported grid type
    if (grid.type() == "unstructured") {
        const atlas::UnstructuredGrid ugrid(grid);
        for (size_t n = 0; n < nb_nodes; ++n) {
            if (ghost[n] == 0) {
                ugrid.lonlat(global_index[n] - global_index_base, lonlat);
                if (!equal_within_roundoff(lonlat[0], lon[n]) || !equal_within_roundoff(lonlat[1], lat[n])) {
                    throw_Exception("In MeshBuilder: UnstructuredGrid from config does not match mesh coordinates", Here());
                }
            }
        }
    } else if (grid.type() == "structured") {
        const atlas::StructuredGrid sgrid(grid);
        for (size_t n = 0; n < nb_nodes; ++n) {
            if (ghost[n] == 0) {
                idx_t i, j;
                sgrid.index2ij(global_index[n] - global_index_base, i, j);
                sgrid.lonlat(i, j, lonlat);
                if (!equal_within_roundoff(lonlat[0], lon[n]) || !equal_within_roundoff(lonlat[1], lat[n])) {
                    throw_Exception("In MeshBuilder: StructuredGrid from config does not match mesh coordinates", Here());
                }
            }
        }
    } else {
        for (size_t n = 0; n < nb_nodes; ++n) {
            if (ghost[n] == 0) {
                auto point = *(grid.lonlat().begin() + static_cast<size_t>(global_index[n] - global_index_base));
                if (!equal_within_roundoff(point.lon(), lon[n]) || !equal_within_roundoff(point.lat(), lat[n])) {
                    throw_Exception("In MeshBuilder: Grid from config does not match mesh coordinates", Here());
                }
            }
        }
    }
}
}  // namespace detail

namespace mesh {

//----------------------------------------------------------------------------------------------------------------------

Mesh MeshBuilder::operator()(const std::vector<double>& lon, const std::vector<double>& lat,
                             const std::vector<int>& ghost, const std::vector<gidx_t>& global_index,
                             const std::vector<idx_t>& remote_index, const idx_t remote_index_base,
                             const std::vector<int>& partition,
                             const std::vector<std::array<gidx_t, 3>>& triag_nodes_global,
                             const std::vector<gidx_t>& triag_global_index,
                             const std::vector<std::array<gidx_t, 4>>& quad_nodes_global,
                             const std::vector<gidx_t>& quad_global_index,
                             const eckit::Configuration& config) const {
    constexpr gidx_t global_index_base = 1;

    return operator()(
        span<const gidx_t>{global_index.data(),global_index.size()},
        span<const double>{lon.data(), lon.size()},  span<const double>{lat.data(), lat.size()},
        span<const double>{lon.data(), lon.size()},  span<const double>{lat.data(), lat.size()},
        span<const int>{ghost.data(), ghost.size()}, span<const int>{partition.data(), partition.size()},
        span<const idx_t>{remote_index.data(), remote_index.size()}, remote_index_base,
        span<const gidx_t>{triag_global_index.data(), triag_global_index.size()},
        mdspan<const gidx_t,extents<size_t,dynamic_extent,3>>{reinterpret_cast<const gidx_t*>(triag_nodes_global.data()),triag_nodes_global.size()},
        span<const gidx_t>{quad_global_index.data(), quad_global_index.size()},
        mdspan<const gidx_t,extents<size_t,dynamic_extent,4>>{reinterpret_cast<const gidx_t*>(quad_nodes_global.data()),quad_nodes_global.size()},
        global_index_base,
        config);
}

Mesh MeshBuilder::operator()(size_t nb_nodes, const double lon[], const double lat[], const int ghost[],
                             const gidx_t global_index[], const idx_t remote_index[], const idx_t remote_index_base,
                             const int partition[], size_t nb_triags, const gidx_t triag_nodes_global[],
                             const gidx_t triag_global_index[], size_t nb_quads, const gidx_t quad_nodes_global[],
                             const gidx_t quad_global_index[],
                             const eckit::Configuration& config) const {
    constexpr gidx_t global_index_base = 1;
    return this->operator()(
        nb_nodes, global_index,
        lon, lat,
        ghost, partition, remote_index, remote_index_base,
        nb_triags, triag_global_index, triag_nodes_global,
        nb_quads, quad_global_index, quad_nodes_global,
        global_index_base,
        config
    );
}

Mesh MeshBuilder::operator()(
    size_t nb_nodes, const gidx_t global_index[],
    const double lon[], const double lat[],
    const int ghost[], const int partition[], const idx_t remote_index[], const idx_t remote_index_base,
    size_t nb_triags, const gidx_t triag_nodes_global[], const gidx_t triag_global_index[],
    size_t nb_quads,  const gidx_t quad_nodes_global[],  const gidx_t quad_global_index[],
    gidx_t global_index_base,
    const eckit::Configuration& config) const {
    return this->operator()(
        nb_nodes, global_index,
        lon, lat, 1, 1,
        lon, lat, 1, 1,
        ghost, partition, remote_index, remote_index_base,
        nb_triags, triag_global_index, triag_nodes_global,
        nb_quads, quad_global_index, quad_nodes_global,
        global_index_base,
        config
    );
}


Mesh MeshBuilder::operator()(size_t nb_nodes, const gidx_t global_index[],
                             const double x[], const double y[], size_t xstride, size_t ystride,
                             const double lon[], const double lat[], size_t lonstride, size_t latstride,
                             const int ghost[], const int partitions[], const idx_t remote_index[], const idx_t remote_index_base,
                             size_t nb_triags, const gidx_t triag_global_index[], const gidx_t triag_nodes_global[],
                             size_t nb_quads,  const gidx_t quad_global_index[],  const gidx_t quad_nodes_global[],
                             const eckit::Configuration& config) const {
    constexpr gidx_t global_index_base = 1;
    return this->operator()(
        nb_nodes, global_index,
        x, y, xstride, ystride,
        lon, lat, lonstride, latstride,
        ghost, partitions, remote_index, remote_index_base,
        nb_triags, triag_global_index, triag_nodes_global,
        nb_quads, quad_global_index, quad_nodes_global,
        global_index_base,
        config);
}


Mesh MeshBuilder::operator()(size_t nb_nodes, const gidx_t global_index[],
                             const double x[], const double y[], size_t xstride, size_t ystride,
                             const double lon[], const double lat[], size_t lonstride, size_t latstride,
                             const int ghost[], const int partitions[], const idx_t remote_index[], const idx_t remote_index_base,
                             size_t nb_triags, const gidx_t triag_global_index[], const gidx_t triag_nodes_global[],
                             size_t nb_quads,  const gidx_t quad_global_index[],  const gidx_t quad_nodes_global[],
                             gidx_t global_index_base,
                             const eckit::Configuration& config) const {
    // Get MPI comm from config name or fall back to atlas default comm
    auto mpi_comm_name = [](const auto& config) {
        return config.getString("mpi_comm", atlas::mpi::comm().name());
    };

    Mesh mesh{};
    mesh.metadata().set("mpi_comm", mpi_comm_name(config));
    auto& comm = mpi::comm(mesh.mpi_comm());

    // Setup a grid, if requested via config argument
    if (config.has("grid")) {
        atlas::Grid grid;
        if (config.has("grid.type") && (config.getString("grid.type") == "unstructured")
            && !config.has("grid.xy")) {
            // Assemble the unstructured grid by gathering input lons,lats across ranks
            grid = ::atlas::detail::assemble_unstructured_grid(nb_nodes, lon, lat, ghost, comm);
        } else {
            // Build grid directly from config
            grid = atlas::Grid(config.getSubConfiguration("grid"));
            const bool validate = config.getBool("validate", false);
            if (validate) {
                ::atlas::detail::validate_grid_vs_mesh(comm, grid, nb_nodes, lon, lat, lonstride, latstride, ghost, global_index, global_index_base);
            }
        }
        mesh.setGrid(grid);
    }

    // Populate node data

    mesh.nodes().resize(nb_nodes);
    auto xy        = array::make_view<double, 2>(mesh.nodes().xy());
    auto lonlat    = array::make_view<double, 2>(mesh.nodes().lonlat());
    auto ghostv    = array::make_view<int, 1>(mesh.nodes().ghost());
    auto gidx      = array::make_view<gidx_t, 1>(mesh.nodes().global_index());
    auto ridx      = array::make_indexview<idx_t, 1>(mesh.nodes().remote_index());
    auto partition = array::make_view<int, 1>(mesh.nodes().partition());
    auto halo      = array::make_view<int, 1>(mesh.nodes().halo());

    for (size_t i = 0; i < nb_nodes; ++i) {
        xy(i, size_t(XX))      = x[i*xstride];
        xy(i, size_t(YY))      = y[i*ystride];
        lonlat(i, size_t(LON)) = lon[i*lonstride];
        lonlat(i, size_t(LAT)) = lat[i*latstride];
        ghostv(i)              = ghost[i];
        gidx(i)                = global_index[i] - global_index_base + 1; // Make 1-based!
        ridx(i)                = remote_index[i] - remote_index_base;
        partition(i)           = partitions[i];
    }
    halo.assign(0);

    // Populate cell/element data

    // First, count how many cells of each type are on this processor
    // Then optimize away the element type if globally nb_triags or nb_quads is zero
    size_t sum_nb_triags = 0;
    comm.allReduce(nb_triags, sum_nb_triags, eckit::mpi::sum());
    const bool add_triags = (sum_nb_triags > 0);

    size_t sum_nb_quads = 0;
    comm.allReduce(nb_quads, sum_nb_quads, eckit::mpi::sum());
    const bool add_quads = (sum_nb_quads > 0);

    if (add_triags) {
        mesh.cells().add(mesh::ElementType::create("Triangle"), nb_triags);
    }
    if (add_quads) {
        mesh.cells().add(mesh::ElementType::create("Quadrilateral"), nb_quads);
    }

    atlas::mesh::HybridElements::Connectivity& node_connectivity = mesh.cells().node_connectivity();
    auto cells_part = array::make_view<int, 1>(mesh.cells().partition());
    auto cells_gidx = array::make_view<gidx_t, 1>(mesh.cells().global_index());

    std::unordered_map<gidx_t,idx_t> global_to_local_index;
    global_to_local_index.reserve(nb_nodes);
    for (gidx_t j=0; j<nb_nodes; ++j) {
        global_to_local_index.insert({global_index[j],j});
    }

    size_t idx = 0;
    if (add_triags) {
        idx_t buffer[3];
        for (size_t tri = 0; tri < nb_triags; ++tri) {
            for (size_t i = 0; i < 3; ++i) {
                buffer[i] = global_to_local_index.at(triag_nodes_global[3 * tri + i]);
            }
            node_connectivity.set(idx, buffer);
            cells_gidx(idx) = triag_global_index[tri] - global_index_base + 1; // make 1-based
            idx++;
        }
    }
    if (add_quads) {
        idx_t buffer[4];
        for (size_t quad = 0; quad < nb_quads; ++quad) {
            for (size_t i = 0; i < 4; ++i) {
                buffer[i] = global_to_local_index.at(quad_nodes_global[4 * quad + i]);
            }
            node_connectivity.set(idx, buffer);
            cells_gidx(idx) = quad_global_index[quad] - global_index_base + 1; // make 1-based
            idx++;
        }
    }

    ATLAS_ASSERT(idx == nb_triags + nb_quads);

    cells_part.assign(comm.rank());

    return mesh;
}

Mesh MeshBuilder::operator()(
    span<const gidx_t> global_index,
    strided_span<const double> x,   strided_span<const double> y,
    strided_span<const double> lon, strided_span<const double> lat,
    span<const int> ghost, span<const int> partition,
    span<const idx_t> remote_index, const idx_t remote_index_base,
    span<const gidx_t> triag_global_index, mdspan<const gidx_t, extents<size_t,dynamic_extent,3>> triag_nodes_global,
    span<const gidx_t> quad_global_index,  mdspan<const gidx_t, extents<size_t,dynamic_extent,4>> quad_nodes_global,
    const gidx_t global_index_base,
    const eckit::Configuration& config) const {

    const size_t nb_nodes  = global_index.size();
    const size_t nb_triags = triag_global_index.size();
    const size_t nb_quads  = quad_global_index.size();

    ATLAS_ASSERT(nb_nodes  == lon.size());
    ATLAS_ASSERT(nb_nodes  == lat.size());
    ATLAS_ASSERT(nb_nodes  == ghost.size());
    ATLAS_ASSERT(nb_nodes  == remote_index.size());
    ATLAS_ASSERT(nb_nodes  == partition.size());
    ATLAS_ASSERT(nb_triags == triag_nodes_global.extent(0));
    ATLAS_ASSERT(nb_quads  == quad_nodes_global.extent(0));

    return operator()(
        nb_nodes, global_index.data_handle(),
        x.data_handle(),   y.data_handle(),   x.stride(0),   y.stride(0),
        lon.data_handle(), lat.data_handle(), lon.stride(0), lat.stride(0),
        ghost.data_handle(), partition.data_handle(), remote_index.data_handle(), remote_index_base,
        nb_triags, triag_global_index.data_handle(), triag_nodes_global.data_handle(),
        nb_quads,  quad_global_index.data_handle(),  quad_nodes_global.data_handle(),
        global_index_base,
        config);
}

//----------------------------------------------------------------------------------------------------------------------

Mesh TriangularMeshBuilder::operator()(size_t nb_nodes,  const gidx_t nodes_global_index[], const double x[], const double y[], const double lon[], const double lat[],
                                       size_t nb_triags, const gidx_t triangle_global_index[], const gidx_t triangle_nodes_global_index[]) const {
    constexpr gidx_t global_index_base = 1;
    constexpr size_t stride_1 = 1;
    return this->operator()(
        nb_nodes, nodes_global_index, x, y, stride_1, stride_1, lon, lat, stride_1, stride_1,
        nb_triags, triangle_global_index, triangle_nodes_global_index,
        global_index_base
    );
}

//----------------------------------------------------------------------------------------------------------------------

Mesh TriangularMeshBuilder::operator()(size_t nb_nodes,  const gidx_t nodes_global_index[],
                                       const double x[], const double y[], size_t xstride, size_t ystride,
                                       const double lon[], const double lat[], size_t lonstride, size_t latstride,
                                       size_t nb_triags, const gidx_t triangle_global_index[], const gidx_t triangle_nodes_global_index[],
                                       gidx_t global_index_base) const {
    std::vector<int> ghost(nb_nodes,0);
    std::vector<int> partition(nb_nodes,0);
    std::vector<idx_t> remote_index(nb_nodes);
    idx_t remote_index_base = 0;
    std::iota(remote_index.begin(), remote_index.end(), remote_index_base);

    size_t nb_quads = 0;
    gidx_t quad_global_index[] = {};
    gidx_t quad_nodes_global_index[] = {};

    return meshbuilder_(
        nb_nodes, nodes_global_index,
        x, y, xstride, ystride, lon, lat, lonstride, latstride,
        ghost.data(), partition.data(), remote_index.data(), remote_index_base,
        nb_triags, triangle_global_index, triangle_nodes_global_index,
        nb_quads, quad_global_index, quad_nodes_global_index,
        global_index_base);
}

//----------------------------------------------------------------------------------------------------------------------

}  // namespace mesh
}  // namespace atlas
