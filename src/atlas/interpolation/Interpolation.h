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

#include "eckit/memory/SharedPtr.h"

#include "atlas/interpolation/method/Method.h"

namespace eckit {
class Parametrisation;
}

namespace atlas {
class Field;
class FieldSet;
class FunctionSpace;
class Grid;
namespace interpolation {
class Method;
}
}  // namespace atlas

namespace atlas {

class Interpolation {
public:
    using Implementation = interpolation::Method;
    using Config         = eckit::Parametrisation;

    Interpolation() {}
    Interpolation( const Interpolation& );

    // Setup Interpolation from source to target function space
    Interpolation( const Config&, const FunctionSpace& source, const FunctionSpace& target );

    // Setup Interpolation from source to coordinates given in a field with multiple components
    Interpolation( const Config&, const FunctionSpace& source, const Field& target );

    // Setup Interpolation from source to coordinates given by separate fields for each component
    Interpolation( const Config&, const FunctionSpace& source, const FieldSet& target );

    // Setup Interpolation from source grid to target grid
    Interpolation( const Config&, const Grid& source, const Grid& target );

    void execute( const FieldSet& source, FieldSet& target ) const;

    void execute( const Field& source, Field& target ) const;

    const Implementation* get() const;

    operator bool() const;

    void print( std::ostream& out ) const;

    const FunctionSpace& source() const;
    const FunctionSpace& target() const;

private:
    eckit::SharedPtr<Implementation> implementation_;
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

Interpolation::Implementation* atlas__Interpolation__new_tgt_field( const eckit::Parametrisation* config,
                                                                    const functionspace::FunctionSpaceImpl* source,
                                                                    const field::FieldImpl* target );

Interpolation::Implementation* atlas__Interpolation__new_tgt_fieldset( const eckit::Parametrisation* config,
                                                                       const functionspace::FunctionSpaceImpl* source,
                                                                       const field::FieldSetImpl* target );

void atlas__Interpolation__delete( Interpolation::Implementation* This );
void atlas__Interpolation__execute_field( Interpolation::Implementation* This, const field::FieldImpl* source,
                                          field::FieldImpl* target );
void atlas__Interpolation__execute_fieldset( Interpolation::Implementation* This, const field::FieldSetImpl* source,
                                             field::FieldSetImpl* target );
}

}  // namespace atlas
