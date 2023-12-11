/*
 * (C) Copyright 2023 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#pragma once

#include "atlas/mesh/MeshBuilder.h"

//----------------------------------------------------------------------------------------------------------------------

namespace atlas {
namespace mesh {

// C wrapper interfaces to C++ routines
extern "C" {
TriangularMeshBuilder* atlas__TriangularMeshBuilder__new();
void atlas__TriangularMeshBuilder__delete(TriangularMeshBuilder* This);

Mesh::Implementation* atlas__TriangularMeshBuilder__operator(TriangularMeshBuilder* This,
  size_t nb_nodes, const gidx_t node_global_index[], const double x[], const double y[], const double lon[], const double lat[],
  size_t nb_triags, const gidx_t triangle_global_index[], const gidx_t triangle_nodes_global_index[]);
}

//----------------------------------------------------------------------------------------------------------------------

}  // namespace mesh
}  // namespace atlas
