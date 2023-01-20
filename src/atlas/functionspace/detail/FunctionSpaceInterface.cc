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
#include "atlas/runtime/Exception.h"
#include "atlas/util/Config.h"

namespace atlas {
namespace functionspace {

//-----------------------------------------------------------------------------

// C wrapper interfaces to C++ routines
extern "C" {
void atlas__FunctionSpace__delete(FunctionSpaceImpl* This) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_FunctionSpace");
    delete This;
    This = nullptr;
}

void atlas__FunctionSpace__name(const FunctionSpaceImpl* This, char*& name, int& size) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_FunctionSpace");
    std::string s = This->type();
    size          = static_cast<int>(s.size());
    name          = new char[size + 1];
    std::strncpy(name, s.c_str(), size + 1);
}

field::FieldImpl* atlas__FunctionSpace__create_field(const FunctionSpaceImpl* This,
                                                     const eckit::Configuration* options) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_FunctionSpace");
    ATLAS_ASSERT(options != nullptr);
    field::FieldImpl* field;
    {
        Field f = This->createField(*options);
        field   = f.get();
        field->attach();
    }
    field->detach();
    return field;
}

//------------------------------------------------------------------------------

field::FieldImpl* atlas__FunctionSpace__create_field_template(const FunctionSpaceImpl* This,
                                                              const field::FieldImpl* field_template,
                                                              const eckit::Configuration* options) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_FunctionSpace");
    ATLAS_ASSERT(field_template != nullptr, "Cannot access uninitialised atlas_Field");
    ATLAS_ASSERT(options != nullptr);
    field::FieldImpl* field;
    {
        Field f = This->createField(Field(field_template), *options);
        field   = f.get();
        field->attach();
    }
    field->detach();
    return field;
}

//------------------------------------------------------------------------------

void atlas__FunctionSpace__halo_exchange_field(const FunctionSpaceImpl* This, field::FieldImpl* field) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_FunctionSpace");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    Field f(field);
    This->haloExchange(f);
}

//------------------------------------------------------------------------------

void atlas__FunctionSpace__halo_exchange_fieldset(const FunctionSpaceImpl* This,
                                                  field::FieldSetImpl* fieldset) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_FunctionSpace");
    ATLAS_ASSERT(fieldset != nullptr, "Cannot access uninitialised atlas_FieldSet");
    FieldSet f(fieldset);
    This->haloExchange(f);
}

//------------------------------------------------------------------------------

void atlas__FunctionSpace__adjoint_halo_exchange_field(const FunctionSpaceImpl* This,
                                                       field::FieldImpl* field) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_FunctionSpace");
    ATLAS_ASSERT(field != nullptr, "Cannot access uninitialised atlas_Field");
    Field f(field);
    This->adjointHaloExchange(f);
}

//------------------------------------------------------------------------------

void atlas__FunctionSpace__adjoint_halo_exchange_fieldset(const FunctionSpaceImpl* This,
                                                          field::FieldSetImpl* fieldset) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_FunctionSpace");
    ATLAS_ASSERT(fieldset != nullptr, "Cannot access uninitialised atlas_FieldSet");
    FieldSet f(fieldset);
    This->adjointHaloExchange(f);
}
}

// ------------------------------------------------------------------

}  // namespace functionspace

// ------------------------------------------------------------------

}  // namespace atlas
