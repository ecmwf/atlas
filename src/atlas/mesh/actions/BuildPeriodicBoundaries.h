/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef BuildPeriodicBoundaries_h
#define BuildPeriodicBoundaries_h

namespace atlas {
class Mesh;

namespace mesh {
namespace actions {
/*
 * Make the mesh periodic
 * - Find out which nodes at the WEST-BC are master nodes
 *   of nodes at the EAST-BC
 * - The remote_idx and partition of the EAST-BC nodes
 *   are set to the master nodes at WEST-BC, whereas the
 *   global index remains unchanged
 */
void build_periodic_boundaries( Mesh& mesh );

// ------------------------------------------------------------------
// C wrapper interfaces to C++ routines
extern "C"
{
  void atlas__build_periodic_boundaries (Mesh::mesh_t* mesh);
}
// ------------------------------------------------------------------

} // namespace actions
} // namespace mesh
} // namespace atlas

#endif // BuildPeriodicBoundaries_h
