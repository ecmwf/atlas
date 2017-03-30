/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef BuildEdges_h
#define BuildEdges_h

#include <string>

namespace atlas {
class Mesh;

namespace mesh {
namespace actions {

void build_edges( Mesh& mesh );
void build_pole_edges( Mesh& mesh );
void build_element_to_edge_connectivity( Mesh& mesh );
void build_node_to_edge_connectivity( Mesh& mesh );

// ------------------------------------------------------------------
// C wrapper interfaces to C++ routines
extern "C"
{
  void atlas__build_edges (Mesh::mesh_t* mesh);
  void atlas__build_pole_edges (Mesh::mesh_t* mesh);
  void atlas__build_node_to_edge_connectivity (Mesh::mesh_t* mesh);
}
// ------------------------------------------------------------------

} // namespace actions
} // namespace mesh
} // namespace atlas

#endif // BuildEdges_h
