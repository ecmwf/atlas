/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/mesh/ElementType.h"
#include "atlas/runtime/ErrorHandling.h"
#include "atlas/util/CoordinateEnums.h"

namespace atlas {
namespace mesh {

//------------------------------------------------------------------------------

ElementType* ElementType::create( const std::string& ) {
    NOTIMP;
}

ElementType::ElementType() {}
ElementType::~ElementType() {}

//-----------------------------------------------------------------------------

extern "C" {
void atlas__mesh__ElementType__delete( ElementType* This ) {
    ATLAS_ERROR_HANDLING( delete This );
}
ElementType* atlas__mesh__Triangle__create() {
    return new temporary::Triangle();
}
ElementType* atlas__mesh__Quadrilateral__create() {
    return new temporary::Quadrilateral();
}
ElementType* atlas__mesh__Line__create() {
    return new temporary::Line();
}

const char* atlas__mesh__ElementType__name( const ElementType* This ) {
    ATLAS_ERROR_HANDLING( ASSERT( This ); return This->name().c_str(); );
    return 0;
}

size_t atlas__mesh__ElementType__nb_nodes( const ElementType* This ) {
    ATLAS_ERROR_HANDLING( ASSERT( This ); return This->nb_nodes() );
    return 0;
}

size_t atlas__mesh__ElementType__nb_edges( const ElementType* This ) {
    ATLAS_ERROR_HANDLING( ASSERT( This ); return This->nb_edges() );
    return 0;
}

int atlas__mesh__ElementType__parametric( const ElementType* This ) {
    ATLAS_ERROR_HANDLING( ASSERT( This ); return This->parametric() );
    return 0;
}
}

}  // namespace mesh
}  // namespace atlas
