/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

/// @file   BuildParallelFields.h
/// @author Willem Deconinck
/// @date   June 2014

#pragma once

#include "atlas/library/config.h"
#include "atlas/mesh/Mesh.h"

namespace atlas {
namespace mesh {
class Nodes;
}
}  // namespace atlas

namespace atlas {
namespace mesh {
namespace actions {

/*
 * Build all parallel fields in the mesh
 *  - calls build_nodes_parallel_fields()
 */
void build_parallel_fields(Mesh& mesh);

/*
 * Build parallel fields for the "nodes" function space if they don't exist.
 * - glb_idx:    create unique indices for non-positive values
 * - partition:  set to mpi::rank() for negative values
 * - remote_idx: rebuild from scratch
 */
void build_nodes_parallel_fields(Mesh& mesh);
void build_nodes_parallel_fields(mesh::Nodes& nodes); // deprecated (WARNING: does not change MPI scope)

/*
 * Build parallel fields for the "edges" function space if they don't exist.
 * - glb_idx:    create unique indices for non-positive values
 * - partition:  set to partition of node with lowest glb_idx
 * - remote_idx: rebuild from scratch
 *
 *  TODO: Make sure that the edge-partition is at least one of the partition
 * numbers of the
 *        neighbouring elements.
 *        Because of this problem, the size of the halo should be set to 2
 * instead of 1!!!
 */
void build_edges_parallel_fields(Mesh& mesh);

void build_cells_parallel_fields(Mesh& mesh);

void renumber_nodes_glb_idx(Mesh& mesh);
void renumber_nodes_glb_idx(mesh::Nodes& nodes); // deprecated (WARNING: does not change MPI scope)

// ------------------------------------------------------------------
// C wrapper interfaces to C++ routines

extern "C" {
void atlas__build_parallel_fields(Mesh::Implementation* mesh);
void atlas__build_nodes_parallel_fields(mesh::Nodes* nodes);
void atlas__build_edges_parallel_fields(Mesh::Implementation* mesh);
void atlas__renumber_nodes_glb_idx(mesh::Nodes* nodes);
}

// ------------------------------------------------------------------

}  // namespace actions
}  // namespace mesh
}  // namespace atlas
