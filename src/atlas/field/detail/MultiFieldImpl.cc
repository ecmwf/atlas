/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

// file deepcode ignore CppMemoryLeak: static pointers for global registry are OK and will be cleaned up at end

#include "atlas/field/MultiField.h"
#include "atlas/field/MultiFieldCreator.h"

#include <map>
#include <sstream>
#include <string>

#include "atlas/field/detail/MultiFieldImpl.h"

namespace atlas {
namespace field {

void MultiFieldImpl::add(Field& field) {
    ATLAS_ASSERT(not fieldset_.has(field.name()), "Field with name \"" + field.name() + "\" already exists!");
    fieldset_.add(field);
    MultiFieldArrayRegistry::instance().add(field, array_);
}

}
}
