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

#include "atlas/interpolation/method/Method.h"
#include "atlas/library/config.h"
#include "atlas/util/ObjectHandle.h"

#include "atlas/interpolation/Cache.h"

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

class Interpolation : DOXYGEN_HIDE(public util::ObjectHandle<interpolation::Method>) {
public:
    using Config   = eckit::Parametrisation;
    using Cache    = interpolation::Cache;
    using Metadata = interpolation::Method::Metadata;

    using Handle::Handle;
    Interpolation() = default;

    // Setup Interpolation from source to target function space
    Interpolation(const Config&, const FunctionSpace& source, const FunctionSpace& target) noexcept(false);

    // Setup Interpolation from source to coordinates given in a field with multiple components
    Interpolation(const Config&, const FunctionSpace& source, const Field& target) noexcept(false);

    // Setup Interpolation from source to coordinates given by separate fields for each component
    Interpolation(const Config&, const FunctionSpace& source, const FieldSet& target) noexcept(false);

    // Setup Interpolation from source grid to target grid
    Interpolation(const Config&, const Grid& source, const Grid& target) noexcept(false);

    Metadata execute(const FieldSet& source, FieldSet& target) const;

    Metadata execute(const Field& source, Field& target) const;

    Metadata execute_adjoint(FieldSet& source, const FieldSet& target) const;

    Metadata execute_adjoint(Field& source, const Field& target) const;

    void print(std::ostream& out) const;

    const FunctionSpace& source() const;
    const FunctionSpace& target() const;

    Cache createCache() const;

    Interpolation(const Config&, const Grid& source, const Grid& target, const Cache&) noexcept(false);

    Interpolation(const Config&, const FunctionSpace& source, const FunctionSpace& target, const Cache&) noexcept(false);

    friend std::ostream& operator<<(std::ostream& out, const Interpolation& i) {
        i.print(out);
        return out;
    }
};

#ifndef DOXYGEN_SHOULD_SKIP_THIS
namespace functionspace {
class FunctionSpaceImpl;
}
namespace field {
class FieldImpl;
class FieldSetImpl;
}  // namespace field

extern "C" {

Interpolation::Implementation* atlas__Interpolation__new(const eckit::Parametrisation* config,
                                                         const functionspace::FunctionSpaceImpl* source,
                                                         const functionspace::FunctionSpaceImpl* target);

Interpolation::Implementation* atlas__Interpolation__new_tgt_field(const eckit::Parametrisation* config,
                                                                   const functionspace::FunctionSpaceImpl* source,
                                                                   const field::FieldImpl* target);

Interpolation::Implementation* atlas__Interpolation__new_tgt_fieldset(const eckit::Parametrisation* config,
                                                                      const functionspace::FunctionSpaceImpl* source,
                                                                      const field::FieldSetImpl* target);

void atlas__Interpolation__delete(Interpolation::Implementation* This);
void atlas__Interpolation__execute_field(Interpolation::Implementation* This, const field::FieldImpl* source,
                                         field::FieldImpl* target);
void atlas__Interpolation__execute_fieldset(Interpolation::Implementation* This, const field::FieldSetImpl* source,
                                            field::FieldSetImpl* target);
void atlas__Interpolation__execute_adjoint_field(Interpolation::Implementation* This, field::FieldImpl* source,
                                                 const field::FieldImpl* target);
void atlas__Interpolation__execute_adjoint_fieldset(Interpolation::Implementation* This, field::FieldSetImpl* source,
                                                    const field::FieldSetImpl* target);
}
#endif

}  // namespace atlas
