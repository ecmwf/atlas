/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor
 * does it submit to any jurisdiction.
 */

#include "atlas/mesh/detail/MeshIntf.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/runtime/ErrorHandling.h"

namespace atlas {
namespace mesh {

//----------------------------------------------------------------------------------------------------------------------
// C wrapper interfaces to C++ routines

Mesh::Implementation* atlas__Mesh__new() {
    return new Mesh::Implementation();
}

void atlas__Mesh__delete( Mesh::Implementation* This ) {
    delete This;
}

Nodes* atlas__Mesh__nodes( Mesh::Implementation* This ) {
    ATLAS_ERROR_HANDLING( ASSERT( This ); return &This->nodes(); );
    return nullptr;
}

Edges* atlas__Mesh__edges( Mesh::Implementation* This ) {
    ATLAS_ERROR_HANDLING( ASSERT( This ); return &This->edges(); );
    return nullptr;
}

Cells* atlas__Mesh__cells( Mesh::Implementation* This ) {
    ATLAS_ERROR_HANDLING( ASSERT( This ); return &This->cells(); );
    return nullptr;
}

size_t atlas__Mesh__footprint( Mesh::Implementation* This ) {
    size_t size( 0 );
    ATLAS_ERROR_HANDLING( ASSERT( This ); size = This->footprint(); );
    return size;
}

void atlas__Mesh__clone_to_device( Mesh::Implementation* This ) {
    This->cloneToDevice();
}

void atlas__Mesh__clone_from_device( Mesh::Implementation* This ) {
    This->cloneFromDevice();
}

void atlas__Mesh__sync_host_device( Mesh::Implementation* This ) {
    This->syncHostDevice();
}

//----------------------------------------------------------------------------------------------------------------------

}  // namespace mesh
}  // namespace atlas
