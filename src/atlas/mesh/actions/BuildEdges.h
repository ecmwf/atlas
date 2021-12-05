/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#pragma once

namespace eckit {
class Configuration;
}  // namespace eckit
namespace atlas {
class Mesh;
}  // namespace atlas

namespace atlas {
namespace mesh {
namespace actions {

void build_edges(Mesh& mesh);
void build_edges(Mesh& mesh, const eckit::Configuration&);
void build_pole_edges(Mesh& mesh);
void build_element_to_edge_connectivity(Mesh& mesh);
void build_node_to_edge_connectivity(Mesh& mesh);

// ------------------------------------------------------------------

}  // namespace actions
}  // namespace mesh
}  // namespace atlas

// ------------------------------------------------------------------

namespace atlas {
namespace mesh {
namespace detail {
class MeshImpl;
}  // namespace detail
}  // namespace mesh
}  // namespace atlas


// C wrapper interfaces to C++ routines
extern "C" {
void atlas__build_edges(atlas::mesh::detail::MeshImpl* mesh);
void atlas__build_pole_edges(atlas::mesh::detail::MeshImpl* mesh);
void atlas__build_node_to_edge_connectivity(atlas::mesh::detail::MeshImpl* mesh);
}

// ------------------------------------------------------------------
