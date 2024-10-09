/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

/// @file MultiFieldInterface.h
/// @author Willem Deconinck
/// @date Sep 2014

#pragma once

#include "atlas/field/FieldSet.h"
#include "atlas/field/MultiField.h"

namespace atlas {
namespace functionspace {
}
}  // namespace atlas

namespace atlas {
namespace field {

//----------------------------------------------------------------------------------------------------------------------

// C wrapper interfaces to C++ routines
extern "C" {
MultiFieldImpl* atlas__MultiField__create(eckit::Configuration* config);
MultiFieldImpl* atlas__MultiField__create_shape(int kind, int rank, int shapef[], const char* var_names,
        size_t length, size_t size);
void atlas__MultiField__delete(MultiFieldImpl* This);
int atlas__MultiField__size(MultiFieldImpl* This);
FieldSetImpl* atlas__MultiField__fieldset(MultiFieldImpl* This);
}

//----------------------------------------------------------------------------------------------------------------------

}  // namespace field
}  // namespace atlas
