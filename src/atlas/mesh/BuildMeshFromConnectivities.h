/*
 * (C) Copyright 2023 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include "atlas/mesh/Mesh.h"

#include <array>
#include <vector>

namespace atlas {

//-----------------------------------------------------------------------------

/**
 * \brief Construct a Mesh by importing external connectivity data
 *
 * Given a list of nodes and corresponding cell-to-node connectivity, sets up a Mesh. Not all Mesh
 * fields are initialized, but enough are to build halos and construct a NodeColumns FunctionSpace.
 *
 * The inputs lons, lats, ghosts, global_indices, remote_indices, and partitions are vectors ranging
 * over the nodes locally owned by (or in the ghost nodes of) the MPI task. The global index is a
 * uniform labeling of the nodes across all MPI tasks; the remote index is the 0-based vector index
 * for the node on its owning task.
 *
 * The tri/quad connectivities (boundary_nodes and global_indices) are vectors ranging over the
 * cells owned by the MPI task. Each cell is defined by a list of nodes defining its boundary;
 * note that each boundary node must be locally known (whether as an owned of ghost node on the
 * MPI task), in other words, must be an element of the node global_indices. The cell global index
 * is, here also, a uniform labeling over the of the cells across all MPI tasks.
 *
 * Some limitations of the current implementation (but not inherent):
 * - can only set up triangle and quad cells.
 * - cannot import halos, i.e., cells owned by other MPI tasks; halos can still be subsequently
 *   computed by calling the BuildMesh action.
 * - cannot import node-to-cell connectivity information.
 */
Mesh build_mesh_from_connectivities(const std::vector<double>& lons, const std::vector<double>& lats,
                                    const std::vector<int>& ghosts, const std::vector<gidx_t>& global_indices,
                                    const std::vector<idx_t>& remote_indices, const idx_t remote_index_base,
                                    const std::vector<int>& partitions,
                                    const std::vector<std::array<gidx_t, 3>>& tri_boundary_nodes,
                                    const std::vector<gidx_t>& tri_global_indices,
                                    const std::vector<std::array<gidx_t, 4>>& quad_boundary_nodes,
                                    const std::vector<gidx_t>& quad_global_indices);

//-----------------------------------------------------------------------------

}  // namespace atlas
