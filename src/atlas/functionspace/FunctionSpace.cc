/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/functionspace/FunctionSpace.h"

#include <cassert>
#include <iostream>
#include <limits>
#include <sstream>

#include "eckit/exception/Exceptions.h"
#include "eckit/types/Types.h"

#include "atlas/array/DataType.h"
#include "atlas/field/Field.h"
#include "atlas/field/detail/FieldImpl.h"
#include "atlas/library/config.h"
#include "atlas/mesh/actions/BuildParallelFields.h"
#include "atlas/runtime/ErrorHandling.h"

namespace atlas {
namespace functionspace {

//-----------------------------------------------------------------------------

// C wrapper interfaces to C++ routines
extern "C" {
void atlas__FunctionSpace__delete( FunctionSpaceImpl* This ) {
    ATLAS_ERROR_HANDLING( ASSERT( This ); delete This; This = 0; );
}

void atlas__FunctionSpace__name( const FunctionSpaceImpl* This, char*& name, int& size ) {
    ATLAS_ERROR_HANDLING( ASSERT( This ); std::string s = This->type(); size = s.size() + 1; name = new char[size];
                          strcpy( name, s.c_str() ); );
}

field::FieldImpl* atlas__FunctionSpace__create_field( const FunctionSpaceImpl* This,
                                                      const eckit::Configuration* options ) {
    ATLAS_ERROR_HANDLING( ASSERT( This ); ASSERT( options ); field::FieldImpl * field; {
        Field f = This->createField( *options );
        field   = f.get();
        field->attach();
    } field->detach();
                          return field );
    return 0;
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
}

// ------------------------------------------------------------------

atlas::Field FunctionSpaceImpl::createField( const atlas::Field& field ) const {
    return createField( field, util::NoConfig() );
}

void FunctionSpaceImpl::haloExchange( FieldSet&, bool ) const {
    NOTIMP;
}

void FunctionSpaceImpl::haloExchange( Field&, bool ) const {
    NOTIMP;
}

Field NoFunctionSpace::createField( const eckit::Configuration& ) const {
    NOTIMP;
}
Field NoFunctionSpace::createField( const Field&, const eckit::Configuration& ) const {
    NOTIMP;
}

// ------------------------------------------------------------------

}  // namespace functionspace

// ------------------------------------------------------------------

FunctionSpace::FunctionSpace() : functionspace_( new functionspace::NoFunctionSpace() ) {}

FunctionSpace::FunctionSpace( const Implementation* functionspace ) : functionspace_( functionspace ) {}

FunctionSpace::FunctionSpace( const FunctionSpace& functionspace ) : functionspace_( functionspace.functionspace_ ) {}

std::string FunctionSpace::type() const {
    return functionspace_->type();
}

FunctionSpace::operator bool() const {
    return functionspace_->operator bool();
}

size_t FunctionSpace::footprint() const {
    return functionspace_->footprint();
}

Field FunctionSpace::createField( const eckit::Configuration& config ) const {
    return functionspace_->createField( config );
}

Field FunctionSpace::createField( const Field& other, const eckit::Configuration& config ) const {
    return functionspace_->createField( other, config );
}

std::string FunctionSpace::distribution() const {
    return functionspace_->distribution();
}

void FunctionSpace::haloExchange( Field& field, bool on_device ) const {
    return functionspace_->haloExchange( field, on_device );
}

void FunctionSpace::haloExchange( FieldSet& fields, bool on_device ) const {
    return functionspace_->haloExchange( fields, on_device );
}

// ------------------------------------------------------------------

}  // namespace atlas
