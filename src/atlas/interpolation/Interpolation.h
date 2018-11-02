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

#include "eckit/config/Configuration.h"
#include "eckit/memory/SharedPtr.h"

#include "atlas/interpolation/method/Method.h"

namespace atlas {
class Field;
class FieldSet;
class FunctionSpace;
class Grid;
}  // namespace atlas

namespace atlas {

class Interpolation {
public:
    using Implementation = interpolation::Method;
    using Config         = Implementation::Config;

    Interpolation() {}
    Interpolation( const Interpolation& );
    Interpolation( const Config&, const FunctionSpace& source, const FunctionSpace& target );
    Interpolation( const Config&, const Grid& source, const Grid& target );

    void execute( const FieldSet& source, FieldSet& target ) const { get()->execute( source, target ); }
    void execute( const Field& source, Field& target ) const { get()->execute( source, target ); }

    const Implementation* get() const { return implementation_.get(); }

    operator bool() const { return implementation_; }

    void print( std::ostream& out ) const { implementation_->print( out ); }

    const FunctionSpace& source() const { return implementation_->source(); }
    const FunctionSpace& target() const { return implementation_->target(); }

private:
    eckit::SharedPtr<const Implementation> implementation_;
};

/// C-interface

namespace functionspace {
class FunctionSpaceImpl;
}
namespace field {
class FieldImpl;
class FieldSetImpl;
}  // namespace field

extern "C" {

Interpolation::Implementation* atlas__Interpolation__new( const eckit::Parametrisation* config,
                                                          const functionspace::FunctionSpaceImpl* source,
                                                          const functionspace::FunctionSpaceImpl* target );
void atlas__Interpolation__delete( Interpolation::Implementation* This );
void atlas__Interpolation__execute_field( Interpolation::Implementation* This, const field::FieldImpl* source,
                                          field::FieldImpl* target );
void atlas__Interpolation__execute_fieldset( Interpolation::Implementation* This, const field::FieldSetImpl* source,
                                             field::FieldSetImpl* target );
}

}  // namespace atlas
