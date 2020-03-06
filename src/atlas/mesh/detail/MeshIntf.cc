/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/mesh/detail/MeshIntf.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/runtime/Exception.h"

namespace atlas {
namespace mesh {

//----------------------------------------------------------------------------------------------------------------------
// C wrapper interfaces to C++ routines

Mesh::Implementation* atlas__Mesh__new() {
    return new Mesh::Implementation();
}

void atlas__Mesh__delete( Mesh::Implementation* This ) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialisd atlas_Mesh" );
    delete This;
}

Nodes* atlas__Mesh__nodes( Mesh::Implementation* This ) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialisd atlas_Mesh" );
    return &This->nodes();
}

Edges* atlas__Mesh__edges( Mesh::Implementation* This ) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialisd atlas_Mesh" );
    return &This->edges();
}

Cells* atlas__Mesh__cells( Mesh::Implementation* This ) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialisd atlas_Mesh" );
    return &This->cells();
}

size_t atlas__Mesh__footprint( Mesh::Implementation* This ) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialisd atlas_Mesh" );
    return This->footprint();
}

void atlas__Mesh__update_device( Mesh::Implementation* This ) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialisd atlas_Mesh" );
    This->updateDevice();
}

void atlas__Mesh__update_host( Mesh::Implementation* This ) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialisd atlas_Mesh" );
    This->updateHost();
}

void atlas__Mesh__sync_host_device( Mesh::Implementation* This ) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialisd atlas_Mesh" );
    This->syncHostDevice();
}

//----------------------------------------------------------------------------------------------------------------------

}  // namespace mesh
}  // namespace atlas
