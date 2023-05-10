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
 * over the nodes locally owned by (or in the ghost nodes of) the MPI task. The tri/quad
 * connectivities are vectors ranging over the cells owned by the MPI task. Each cell is defined by
 * a list of nodes defining its boundary; note that each boundary node must be locally known
 * (whether as an owned or ghost node on the MPI task), in other words, must be included in the
 * vectors over nodes.
 *
 * The inputs tri/quad_connectivities give a local and a global index to each cell. The local index
 * is 0-based and counts the cells on the MPI task. The global index is a uniform labeling over the
 * entire grid. The local cell index must range over triangle cells *first* and quad cells *after*.
 *
 * The boundary nodes of each cell are the *local* indices into the *local* vector global_indices.
 * For example, an MPI task may have nodes with global_indices {2,5,18,42,200}, and it may own a
 * triangular cell defined by nodes {2,5,42}, in this case the boundary nodes are the local indices
 * {0,1,3} of the boundary nodes within the vector global_indices.
 *
 * Some restrictions:
 * - currently, only triangle and quad cells are supported (but others could be added)
 * - currently, no support for importing node-to-cell connectivity information (but could be added)
 */
Mesh build_mesh_from_connectivities(const std::vector<double>& lons, const std::vector<double>& lats,
                                    const std::vector<int>& ghosts, const std::vector<gidx_t>& global_indices,
                                    const std::vector<idx_t>& remote_index, const std::vector<int>& partitions,
                                    const std::vector<std::array<gidx_t, 3>>& tri_boundary_nodes,
                                    const std::vector<gidx_t>& tri_global_indices,
                                    const std::vector<std::array<gidx_t, 4>>& quad_boundary_nodes,
                                    const std::vector<gidx_t>& quad_global_indices);

//-----------------------------------------------------------------------------

}  // namespace atlas
