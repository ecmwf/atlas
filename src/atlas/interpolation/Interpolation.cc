/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "eckit/exception/Exceptions.h"

#include "atlas/interpolation/Interpolation.h"
#include "atlas/field/Field.h"
#include "atlas/field/FieldSet.h"
#include "atlas/functionspace/FunctionSpace.h"

namespace atlas {

Interpolation::Interpolation( const Config& config, const FunctionSpace& source, const FunctionSpace& target ) :
    implementation_( [&]() -> Implementation* {
        std::string type;
        config.get( "type", type );
        Implementation* impl = interpolation::MethodFactory::build( type, config );
        impl->setup( source, target );
        return impl;
    }() ) {}

Interpolation::Interpolation( const Interpolation& other ) : implementation_( other.implementation_ ) {}

extern "C" {
Interpolation::Implementation* atlas__Interpolation__new( const eckit::Parametrisation* config,
                                                          const functionspace::FunctionSpaceImpl* source,
                                                          const functionspace::FunctionSpaceImpl* target ) {
    Interpolation::Implementation* interpolator;
    {
        Interpolation im( *config, FunctionSpace( source ), FunctionSpace( target ) );
        interpolator = const_cast<Interpolation::Implementation*>( im.get() );
        interpolator->attach();
    }
    interpolator->detach();
    return interpolator;
}

void atlas__Interpolation__delete( Interpolation::Implementation* This ) {
    delete This;
}

void atlas__Interpolation__execute_field( Interpolation::Implementation* This, const field::FieldImpl* source,
                                          field::FieldImpl* target ) {
    Field t( target );
    This->execute( Field( source ), t );
}

void atlas__Interpolation__execute_fieldset( Interpolation::Implementation* This, const field::FieldSetImpl* source,
                                             field::FieldSetImpl* target ) {
    FieldSet t( target );
    This->execute( FieldSet( source ), t );
}

}  // extern "C"

}  // namespace atlas
