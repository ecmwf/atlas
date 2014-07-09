/*
 * (C) Copyright 1996-2014 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

/// @author Willem Deconinck
/// @date   June 2014

#ifndef BuildParallelFields_hpp
#define BuildParallelFields_hpp
#include <string>
namespace atlas {
  class Mesh;
  class FunctionSpace;
namespace actions {

/*
 * Build all parallel fields in the mesh
 *  - calls build_nodes_parallel_fields()
 */
void build_parallel_fields( Mesh& mesh );

/*
 * Build parallel fields for the "nodes" function space if they don't exist.
 * - glb_idx:    create unique indices for non-positive values
 * - partition:  set to MPL::rank() for negative values
 * - remote_idx: rebuild from scratch
 */
void build_nodes_parallel_fields( FunctionSpace& nodes );

/*
 * Build parallel fields for the "edges" function space if they don't exist.
 * - glb_idx:    create unique indices for non-positive values
 * - partition:  set to partition of node with lowest glb_idx
 * - remote_idx: rebuild from scratch
 *
 *  TODO: Make sure that the edge-partition is at least one of the partition numbers of the
 *        neighbouring elements.
 *        Because of this problem, the size of the halo should be set to 2 instead of 1!!!
 */
void build_edges_parallel_fields( FunctionSpace& edges, FunctionSpace& nodes );

// ------------------------------------------------------------------
// C wrapper interfaces to C++ routines
extern "C"
{
  void atlas__build_parallel_fields (Mesh* mesh);
  void atlas__build_nodes_parallel_fields (FunctionSpace* nodes);
  void atlas__build_edges_parallel_fields (FunctionSpace* edges, FunctionSpace* nodes);
}
// ------------------------------------------------------------------

} // namespace actions
} // namespace atlas


#endif // BuildParallelFields_hpp
