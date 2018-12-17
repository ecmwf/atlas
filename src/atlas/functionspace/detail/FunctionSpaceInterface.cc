/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <cstring>

#include "FunctionSpaceImpl.h"
#include "FunctionSpaceInterface.h"

#include "atlas/field/Field.h"
#include "atlas/field/FieldSet.h"
#include "atlas/field/detail/FieldImpl.h"
#include "atlas/library/config.h"
#include "atlas/runtime/ErrorHandling.h"
#include "atlas/util/Config.h"

namespace atlas {
namespace functionspace {

//-----------------------------------------------------------------------------

// C wrapper interfaces to C++ routines
extern "C" {
void atlas__FunctionSpace__delete( FunctionSpaceImpl* This ) {
    ATLAS_ERROR_HANDLING( ASSERT( This ); delete This; This = nullptr; );
}

void atlas__FunctionSpace__name( const FunctionSpaceImpl* This, char*& name, int& size ) {
    ATLAS_ERROR_HANDLING( ASSERT( This ); std::string s = This->type(); size = static_cast<int>( s.size() + 1 );
                          name = new char[size]; strcpy( name, s.c_str() ); );
}

field::FieldImpl* atlas__FunctionSpace__create_field( const FunctionSpaceImpl* This,
                                                      const eckit::Configuration* options ) {
    ATLAS_ERROR_HANDLING( ASSERT( This ); ASSERT( options ); field::FieldImpl * field; {
        Field f = This->createField( *options );
        field   = f.get();
        field->attach();
    } field->detach();
                          return field );
    return nullptr;
}

//------------------------------------------------------------------------------

field::FieldImpl* atlas__FunctionSpace__create_field_template( const FunctionSpaceImpl* This,
                                                               const field::FieldImpl* field_template,
                                                               const eckit::Configuration* options ) {
    ASSERT( This );
    ASSERT( options );
    field::FieldImpl* field;
    {
        Field f = This->createField( Field( field_template ), *options );
        field   = f.get();
        field->attach();
    }
    field->detach();
    return field;
}

//------------------------------------------------------------------------------

void atlas__FunctionSpace__halo_exchange_field( const FunctionSpaceImpl* This, field::FieldImpl* field ) {
    ATLAS_ERROR_HANDLING( ASSERT( This ); ASSERT( field ); Field f( field ); This->haloExchange( f ); );
}

//------------------------------------------------------------------------------

void atlas__FunctionSpace__halo_exchange_fieldset( const FunctionSpaceImpl* This, field::FieldSetImpl* fieldset ) {
    ATLAS_ERROR_HANDLING( ASSERT( This ); ASSERT( fieldset ); FieldSet f( fieldset ); This->haloExchange( f ); );
}
}

// ------------------------------------------------------------------

}  // namespace functionspace

// ------------------------------------------------------------------

}  // namespace atlas
