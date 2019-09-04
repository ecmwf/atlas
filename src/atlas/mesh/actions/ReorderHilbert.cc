/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/mesh/actions/ReorderHilbert.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"
#include "atlas/runtime/Trace.h"

namespace atlas {
namespace mesh {
namespace actions {

// -------------------------------------------------------------------------------------

ReorderHilbert::ReorderHilbert( Mesh& mesh ) : mesh_( mesh ) {}

void ReorderHilbert::operator()() {
    ATLAS_TRACE( "ReorderHilbert()" );
}

// ------------------------------------------------------------------

}  // namespace actions
}  // namespace mesh
}  // namespace atlas
