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

#include "atlas/mesh/Mesh.h"

//----------------------------------------------------------------------------------------------------------------------

namespace atlas {
namespace mesh {

// C wrapper interfaces to C++ routines
extern "C" {
Mesh::Implementation* atlas__Mesh__new();
void atlas__Mesh__delete(Mesh::Implementation* This);
Nodes* atlas__Mesh__nodes(Mesh::Implementation* This);
Edges* atlas__Mesh__edges(Mesh::Implementation* This);
Cells* atlas__Mesh__cells(Mesh::Implementation* This);
size_t atlas__Mesh__footprint(Mesh::Implementation* This);
void atlas__Mesh__update_device(Mesh::Implementation* This);
void atlas__Mesh__update_host(Mesh::Implementation* This);
void atlas__Mesh__sync_host_device(Mesh::Implementation* This);
}

//----------------------------------------------------------------------------------------------------------------------

}  // namespace mesh
}  // namespace atlas
