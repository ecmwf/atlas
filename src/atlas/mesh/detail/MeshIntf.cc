/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "atlas/mesh/detail/MeshIntf.h"
#include "atlas/runtime/ErrorHandling.h"
#include "atlas/mesh/Nodes.h"

namespace atlas {
namespace mesh {

//----------------------------------------------------------------------------------------------------------------------
// C wrapper interfaces to C++ routines

Mesh::Implementation* atlas__Mesh__new () {
	return new Mesh::Implementation();
}

void atlas__Mesh__delete (Mesh::Implementation* This) {
	delete This;
}

Nodes* atlas__Mesh__create_nodes (Mesh::Implementation* This, int nb_nodes)
{
  ATLAS_ERROR_HANDLING(
    ASSERT( This );
    This->nodes().resize(nb_nodes);
    return &This->nodes();
  );
  return nullptr;
}

Nodes* atlas__Mesh__nodes (Mesh::Implementation* This) {
  ATLAS_ERROR_HANDLING(
    ASSERT( This );
    return &This->nodes();
  );
  return nullptr;
}

Edges* atlas__Mesh__edges (Mesh::Implementation* This) {
  ATLAS_ERROR_HANDLING(
    ASSERT( This );
    return &This->edges();
  );
  return nullptr;
}

Cells* atlas__Mesh__cells (Mesh::Implementation* This) {
  ATLAS_ERROR_HANDLING(
    ASSERT( This );
    return &This->cells();
  );
  return nullptr;
}

size_t atlas__Mesh__footprint (Mesh::Implementation* This) {
  size_t size(0);
  ATLAS_ERROR_HANDLING(
    ASSERT( This );
    size = This->footprint();
  );
  return size;
}

//----------------------------------------------------------------------------------------------------------------------

} // namespace mesh
} // namespace atlas
