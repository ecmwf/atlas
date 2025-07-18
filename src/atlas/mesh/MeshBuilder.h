/*
 * (C) Copyright 2023 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include "eckit/config/Configuration.h"

#include "atlas/parallel/mpi/mpi.h"
#include "atlas/util/Config.h"
#include "atlas/mdspan.h"
#include "atlas/mesh/Mesh.h"

#include <array>
#include <vector>

namespace atlas {
namespace mesh {

//-----------------------------------------------------------------------------

/**
 * \brief Construct a Mesh by importing external connectivity data
 *
 * Given a list of nodes and corresponding cell-to-node connectivity, sets up a Mesh. Not all Mesh
 * fields are initialized, but enough are to build halos and construct a NodeColumns FunctionSpace.
 *
 * Some limitations of the current implementation (but not inherent):
 * - can only set up triangle and quad cells.
 * - cannot import halos, i.e., cells owned by other MPI tasks; halos can still be subsequently
 *   computed by calling the BuildMesh action.
 * - cannot import node-to-cell connectivity information.
 *
 * The Mesh's Grid can be initialized via the call operator's optional config argument.
 */
class MeshBuilder : public eckit::Owned {
public:
    MeshBuilder(const eckit::Configuration& = util::NoConfig()) {}

    /**
     * \brief C-interface to construct a Mesh from external connectivity data
     *
     * The inputs lons, lats, ghosts, global_indices, remote_indices, and partitions are vectors of
     * size nb_nodes, ranging over the nodes locally owned by (or in the ghost nodes of) the MPI
     * task. The global index is a uniform labeling of the nodes across all MPI tasks; the remote
     * index is a remote_index_base-based vector index for the node on its owning task.
     *
     * The tri/quad connectivities (boundary_nodes and global_indices) are vectors ranging over the
     * cells owned by the MPI task. Each cell is defined by a list of nodes defining its boundary;
     * note that each boundary node must be locally known (whether as an owned of ghost node on the
     * MPI task), in other words, must be an element of the node global_indices. The boundary nodes
     * are ordered node-varies-fastest, element-varies-slowest order. The cell global index is,
     * here also, a uniform labeling over the of the cells across all MPI tasks.
     *
     * The config argument can be used to
     * - Request the Mesh's Grid to be constructed, usually from the config. If the Grid is either
     *   an UnstructuredGrid or a Structured grid, the `validate` bool option can be used to trigger
     *   a simple check that the grid is consistent with the lons/lats passed in to the MeshBuilder.
     *   In the special case where `grid.type` is unstructured and the `grid.xy` coordinates are
     *   _not_ given, then the grid is constructed from the lons/lats passed to the MeshBuilder.
     * - Select which MPI communicator to use.
     */

    // operator (1) deprecated, use operator (2)
    [[deprecated("Use operator (2) with extra global_index_base argument")]]
    Mesh operator()(size_t nb_nodes, const double lon[], const double lat[], const int ghost[],
                    const gidx_t global_index[], const idx_t remote_index[], const idx_t remote_index_base,
                    const int partition[],
                    size_t nb_triags, const gidx_t triag_nodes_global[], const gidx_t triag_global_index[],
                    size_t nb_quads,  const gidx_t quad_nodes_global[],  const gidx_t quad_global_index[],
                    const eckit::Configuration& config = util::NoConfig()) const;

    // operator (2) = slightly reordered arguments of (1) and extra global_index_base argument
    Mesh operator()(size_t nb_nodes, const gidx_t global_index[],
                    const double lon[], const double lat[],
                    const int ghost[], const int partition[], const idx_t remote_index[], const idx_t remote_index_base,
                    size_t nb_triags, const gidx_t triag_nodes_global[], const gidx_t triag_global_index[],
                    size_t nb_quads,  const gidx_t quad_nodes_global[],  const gidx_t quad_global_index[],
                    gidx_t global_index_base,
                    const eckit::Configuration& config = util::NoConfig()) const;

    // operator (3) deprecated, use (4)
    // Delegate to next operator with global index base = 1
    [[deprecated("Use operator (4) with extra global_index_base argument")]]
    Mesh operator()(size_t nb_nodes, const gidx_t global_index[],
                    const double x[], const double y[], size_t xstride, size_t ystride,
                    const double lon[], const double lat[], size_t lonstride, size_t latstride,
                    const int ghost[], const int partition[], const idx_t remote_index[], const idx_t remote_index_base,
                    size_t nb_triags, const gidx_t triag_global_index[], const gidx_t triag_nodes_global[],
                    size_t nb_quads,  const gidx_t quad_global_index[],  const gidx_t quad_nodes_global[],
                    const eckit::Configuration& config = util::NoConfig()) const;

    // operator (4) = operator (3) with extra global_index_base argument
    Mesh operator()(size_t nb_nodes, const gidx_t global_index[],
                    const double x[], const double y[], size_t xstride, size_t ystride,
                    const double lon[], const double lat[], size_t lonstride, size_t latstride,
                    const int ghost[], const int partition[], const idx_t remote_index[], const idx_t remote_index_base,
                    size_t nb_triags, const gidx_t triag_global_index[], const gidx_t triag_nodes_global[],
                    size_t nb_quads,  const gidx_t quad_global_index[],  const gidx_t quad_nodes_global[],
                    gidx_t global_index_base,
                    const eckit::Configuration& config = util::NoConfig()) const;

    template<typename T>
    using span = mdspan<T,dims<1>,layout_right>;

    template<typename T>
    using strided_span = mdspan<T,dims<1>,layout_stride>;

    // operator (5)
    Mesh operator()(span<const gidx_t> global_index,
                    strided_span<const double> x,   strided_span<const double> y,
                    strided_span<const double> lon, strided_span<const double> lat,
                    span<const int> ghost, span<const int> partition,
                    span<const idx_t> remote_index, const idx_t remote_index_base,
                    span<const gidx_t> triag_global_index, mdspan<const gidx_t, extents<size_t,dynamic_extent,3>> triag_nodes_global,
                    span<const gidx_t> quad_global_index,  mdspan<const gidx_t, extents<size_t,dynamic_extent,4>> quad_nodes_global,
                    const gidx_t global_index_base,
                    const eckit::Configuration& config = util::NoConfig()) const;

    /**
     * \brief C++-interface to construct a Mesh from external connectivity data
     *
     * Provides a wrapper to the C-interface using STL containers.
     */
    [[deprecated("Use mdspan based operator with extra global_index_base argument")]]
    Mesh operator()(const std::vector<double>& lons, const std::vector<double>& lats, const std::vector<int>& ghosts,
                    const std::vector<gidx_t>& global_index, const std::vector<idx_t>& remote_index,
                    const idx_t remote_index_base, const std::vector<int>& partitions,
                    const std::vector<std::array<gidx_t, 3>>& triag_nodes_global,
                    const std::vector<gidx_t>& triag_global_index,
                    const std::vector<std::array<gidx_t, 4>>& quad_nodes_global,
                    const std::vector<gidx_t>& quad_global_index,
                    const eckit::Configuration& config = util::NoConfig()) const;
};

//-----------------------------------------------------------------------------


class TriangularMeshBuilder {
public:
    TriangularMeshBuilder(const eckit::Configuration& config = util::NoConfig()) :
        meshbuilder_(config) {}

    /**
     * \brief C-interface to construct a Triangular Mesh from external connectivity data
     *
     * The inputs x, y, lons, lats, ghost, global_index, remote_index, and partition are vectors of
     * size nb_nodes, ranging over the nodes locally owned by (or in the ghost nodes of) the MPI
     * task. The global index is a uniform labeling of the nodes across all MPI tasks; the remote
     * index is a remote_index_base-based vector index for the node on its owning task.
     *
     * The triangle connectivities (boundary_nodes and global_index) are vectors ranging over the
     * cells owned by the MPI task. Each cell is defined by a list of nodes defining its boundary;
     * note that each boundary node must be locally known (whether as an owned of ghost node on the
     * MPI task), in other words, must be an element of the node global_indices. The boundary nodes
     * are ordered node-varies-fastest, element-varies-slowest order. The cell global index is,
     * here also, a uniform labeling over the of the cells across all MPI tasks.
     */
    [[deprecated("Use operator with extra global_index_base argument and strides")]]
    Mesh operator()(size_t nb_nodes,    const gidx_t node_global_index[],
                    const double x[],   const double y[],
                    const double lon[], const double lat[],
                    size_t nb_triags,   const gidx_t triangle_global_index[], const gidx_t triangle_nodes_global_index[]) const;

    Mesh operator()(size_t nb_nodes,    const gidx_t node_global_index[],
                    const double x[],   const double y[],   size_t xstride,   size_t ystride,
                    const double lon[], const double lat[], size_t lonstride, size_t latstride,
                    size_t nb_triags,   const gidx_t triangle_global_index[], const gidx_t triangle_nodes_global_index[],
                    gidx_t global_index_base) const;

private:
    MeshBuilder meshbuilder_;
};

//-----------------------------------------------------------------------------

}  // namespace mesh
}  // namespace atlas
