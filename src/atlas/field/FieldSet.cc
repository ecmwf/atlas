/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <sstream>

#include "atlas/field/Field.h"
#include "atlas/field/FieldSet.h"
#include "atlas/grid/Grid.h"
#include "atlas/runtime/Exception.h"

namespace atlas {
namespace field {

//------------------------------------------------------------------------------------------------------

FieldSetImpl::FieldSetImpl(const std::string& /*name*/): name_() {}

void FieldSetImpl::clear() {
    index_.clear();
    fields_.clear();
}

Field FieldSetImpl::add(const Field& field) {
    if (field.name().size()) {
        index_[field.name()] = size();
    }
    else {
        std::stringstream name;
        name << name_ << "[" << size() << "]";
        index_[name.str()] = size();
    }
    fields_.push_back(field);
    return field;
}

bool FieldSetImpl::has(const std::string& name) const {
    return index_.count(name);
}


Field& FieldSetImpl::field(const std::string& name) const {
    if (!has(name)) {
        const std::string msg("FieldSet" + (name_.length() ? " \"" + name_ + "\"" : "") + ": cannot find field \"" +
                              name + "\"");
        throw_Exception(msg, Here());
    }
    return const_cast<Field&>(fields_[index_.at(name)]);
}

void FieldSetImpl::haloExchange(bool on_device) const {
    for (idx_t i = 0; i < size(); ++i) {
        field(i).haloExchange(on_device);
    }
}

void FieldSetImpl::adjointHaloExchange(bool on_device) const {
    for (idx_t i = 0; i < size(); ++i) {
        field(i).adjointHaloExchange(on_device);
    }
}

void FieldSetImpl::set_dirty(bool value) const {
    for (idx_t i = 0; i < size(); ++i) {
        field(i).set_dirty(value);
    }
}

std::vector<std::string> FieldSetImpl::field_names() const {
    std::vector<std::string> ret;

    for (const_iterator field = cbegin(); field != cend(); ++field) {
        ret.push_back(field->name());
    }

    return ret;
}

//-----------------------------------------------------------------------------
// C wrapper interfaces to C++ routines
extern "C" {

FieldSetImpl* atlas__FieldSet__new(char* name) {
    FieldSetImpl* fset = new FieldSetImpl(std::string(name));
    fset->name()       = name;
    return fset;
}

void atlas__FieldSet__delete(FieldSetImpl* This) {
    ATLAS_ASSERT(This != nullptr, "Reason: Use of uninitialised atlas_FieldSet");
    delete This;
}

const char* atlas__FieldSet__name(FieldSetImpl* This) {
    ATLAS_ASSERT(This != nullptr, "Cannot access name of uninitialised atlas_FieldSet");
    return This->name().c_str();
}

void atlas__FieldSet__add_field(FieldSetImpl* This, FieldImpl* field) {
    ATLAS_ASSERT(This != nullptr, "Reason: Use of uninitialised atlas_FieldSet");
    ATLAS_ASSERT(field != nullptr, "Reason: Use of uninitialised atlas_Field");
    This->add(field);
}

int atlas__FieldSet__has_field(const FieldSetImpl* This, char* name) {
    ATLAS_ASSERT(This != nullptr, "Reason: Use of uninitialised atlas_FieldSet");
    return This->has(std::string(name));
}

idx_t atlas__FieldSet__size(const FieldSetImpl* This) {
    ATLAS_ASSERT(This != nullptr, "Reason: Use of uninitialised atlas_FieldSet");
    return This->size();
}

FieldImpl* atlas__FieldSet__field_by_name(FieldSetImpl* This, char* name) {
    ATLAS_ASSERT(This != nullptr, "Reason: Use of uninitialised atlas_FieldSet");
    return This->field(std::string(name)).get();
}

FieldImpl* atlas__FieldSet__field_by_idx(FieldSetImpl* This, idx_t idx) {
    ATLAS_ASSERT(This != nullptr, "Reason: Use of uninitialised atlas_FieldSet");
    return This->operator[](idx).get();
}

void atlas__FieldSet__set_dirty(FieldSetImpl* This, int value) {
    ATLAS_ASSERT(This != nullptr, "Reason: Use of uninitialised atlas_FieldSet");
    This->set_dirty(value);
}

void atlas__FieldSet__halo_exchange(FieldSetImpl* This, int on_device) {
    ATLAS_ASSERT(This != nullptr, "Reason: Use of uninitialised atlas_FieldSet");
    This->haloExchange(on_device);
}
}
//-----------------------------------------------------------------------------

}  // namespace field

//------------------------------------------------------------------------------------------------------

FieldSet::FieldSet(): Handle(new Implementation()) {}

FieldSet::FieldSet(const std::string& name): Handle(new Implementation(name)) {}

FieldSet::FieldSet(const Field& field): Handle(new Implementation()) {
    get()->add(field);
}

const util::Metadata& FieldSet::metadata() const {
    return get()->metadata();
}

util::Metadata& FieldSet::metadata() {
    return get()->metadata();
}

FieldSet FieldSet::clone(const eckit::Parametrisation& config) const {
    FieldSet fset;
    for (idx_t jj = 0; jj < size(); ++jj) {
        fset.add(field(jj).clone(config));
    }
    return fset;
}

void FieldSet::set_dirty(bool value) const {
    get()->set_dirty(value);
}

//------------------------------------------------------------------------------------------------------

}  // namespace atlas
