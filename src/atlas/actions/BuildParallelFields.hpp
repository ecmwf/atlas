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
 * Make the mesh periodic
 * - Find out which nodes at the WEST-BC are master nodes
 *   of nodes at the EAST-BC
 * - The remote_idx and partition of the EAST-BC nodes
 *   are set to the master nodes at WEST-BC, whereas the
 *   global index remains unchanged
 */
void make_periodic( Mesh& mesh );

} // namespace actions
} // namespace atlas

#endif // BuildParallelFields_hpp
