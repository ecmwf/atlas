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


namespace eckit {
class Configuration;
}

namespace atlas {
namespace field {
class FieldImpl;
class FieldSetImpl;
}  // namespace field
}  // namespace atlas

namespace atlas {
namespace functionspace {
class FunctionSpaceImpl;

//------------------------------------------------------------------------------------------------------

// C wrapper interfaces to C++ routines
extern "C" {
void atlas__FunctionSpace__delete(FunctionSpaceImpl* This);
void atlas__FunctionSpace__name(const FunctionSpaceImpl* This, char*& name, int& size);
field::FieldImpl* atlas__FunctionSpace__create_field(const FunctionSpaceImpl* This,
                                                     const eckit::Configuration* options);
field::FieldImpl* atlas__FunctionSpace__create_field_template(const FunctionSpaceImpl* This,
                                                              const field::FieldImpl* field_template,
                                                              const eckit::Configuration* options);
void atlas__FunctionSpace__halo_exchange_field(const FunctionSpaceImpl* This, field::FieldImpl* field);
void atlas__FunctionSpace__halo_exchange_fieldset(const FunctionSpaceImpl* This, field::FieldSetImpl* fieldset);
void atlas__FunctionSpace__adjoint_halo_exchange_field(const FunctionSpaceImpl* This, field::FieldImpl* field);
void atlas__FunctionSpace__adjoint_halo_exchange_fieldset(const FunctionSpaceImpl* This, field::FieldSetImpl* fieldset);
}

//------------------------------------------------------------------------------------------------------

}  // namespace functionspace

//------------------------------------------------------------------------------------------------------

}  // namespace atlas
