/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */


#include <iomanip>
#include <map>
#include <memory>
#include <string>
#include <mutex>

#include "atlas/field/MultiField.h"
#include "atlas/field/detail/MultiFieldImpl.h"
#include "atlas/runtime/Exception.h"

namespace atlas {
namespace field {

//-----------------------------------------------------------------------------

const Field& MultiField::field(const std::string& name) const { return get()->field(name); }
Field& MultiField::field(const std::string& name) { return get()->field(name); }
bool MultiField::has(const std::string& name) const { return get()->has(name); }
std::vector<std::string> MultiField::field_names() const { return get()->field_names(); }

const Field& MultiField::field(const idx_t idx) const { return get()->field(idx); }
Field& MultiField::field(const idx_t idx) { return get()->field(idx); }
idx_t MultiField::size() const { return get()->size(); }

const Field& MultiField::operator[](const idx_t idx) const { return get()->field(idx); }
Field& MultiField::operator[](const idx_t idx) { return get()->field(idx); }

const Field& MultiField::operator[](const std::string& name) const { return get()->field(name); }
Field& MultiField::operator[](const std::string& name) { return get()->field(name); }

const util::Metadata& MultiField::metadata() const { return get()->metadata(); }
util::Metadata& MultiField::metadata() { return get()->metadata(); }

MultiField::operator const array::Array&() const { return get()->array(); }
MultiField::operator array::Array&() { return get()->array(); }

MultiField::operator const FieldSet&() const { return get()->fieldset_; }
MultiField::operator FieldSet&() { return get()->fieldset_; }

const array::Array& MultiField::array() const { return get()->array(); }
array::Array& MultiField::array() { return get()->array(); }

//-----------------------------------------------------------------------------

}  // namespace field
}  // namespace atlas
