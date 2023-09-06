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
#include "atlas/field/detail/FieldInterface.h"
#include "atlas/field/FieldSet.h"
#include "atlas/grid/Grid.h"
#include "atlas/runtime/Exception.h"

namespace atlas {
namespace field {

//------------------------------------------------------------------------------------------------------

void FieldSetImpl::FieldObserver::onFieldRename(FieldImpl& field) {
    std::string name = field.name();
    for (auto& kv: fieldset_.index_) {
        const auto old_name = kv.first;
        const auto idx      = kv.second;
        if (&field == fieldset_.fields_[idx].get()) {
            if (name.empty()) {
                std::stringstream ss;
                ss << fieldset_.name_ << "[" << idx << "]";
                name = ss.str();
            }
            fieldset_.index_.erase(old_name);
            fieldset_.index_[name] = idx;
            return;
        }
    }
    throw_AssertionFailed("Should not be here",Here());
}


FieldSetImpl::FieldSetImpl(const std::string& /*name*/): name_(), field_observer_(*this) {}
FieldSetImpl::~FieldSetImpl() {
    clear();
}

void FieldSetImpl::clear() {
    for( auto& field : fields_ ) {
        field->detachObserver(field_observer_);
    }
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

    field.get()->attachObserver(field_observer_);
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

void atlas__FieldSet__data_int_specf(FieldSetImpl* This, char* name, int*& data, int& rank, int*& shapef, int*& stridesf) {
    atlas__Field__data_int_specf(This->field(name).get(), data, rank, shapef, stridesf);
}

void atlas__FieldSet__data_long_specf(FieldSetImpl* This, char* name, long*& data, int& rank, int*& shapef, int*& stridesf) {
    atlas__Field__data_long_specf(This->field(name).get(), data, rank, shapef, stridesf);
}

void atlas__FieldSet__data_float_specf(FieldSetImpl* This, char* name, float*& data, int& rank, int*& shapef, int*& stridesf) {
    atlas__Field__data_float_specf(This->field(name).get(), data, rank, shapef, stridesf);
}

void atlas__FieldSet__data_double_specf(FieldSetImpl* This, char* name, double*& data, int& rank, int*& shapef, int*& stridesf) {
    atlas__Field__data_double_specf(This->field(name).get(), data, rank, shapef, stridesf);
}

void atlas__FieldSet__data_int_specf_by_idx(FieldSetImpl* This, int& idx, int*& data, int& rank, int*& shapef, int*& stridesf) {
    atlas__Field__data_int_specf(This->operator[](idx).get(), data, rank, shapef, stridesf);
}

void atlas__FieldSet__data_long_specf_by_idx(FieldSetImpl* This, int& idx, long*& data, int& rank, int*& shapef, int*& stridesf) {
    atlas__Field__data_long_specf(This->operator[](idx).get(), data, rank, shapef, stridesf);
}

void atlas__FieldSet__data_float_specf_by_idx(FieldSetImpl* This, int& idx, float*& data, int& rank, int*& shapef, int*& stridesf) {
    atlas__Field__data_float_specf(This->operator[](idx).get(), data, rank, shapef, stridesf);
}

void atlas__FieldSet__data_double_specf_by_idx(FieldSetImpl* This, int& idx, double*& data, int& rank, int*& shapef, int*& stridesf) {
    atlas__Field__data_double_specf(This->operator[](idx).get(), data, rank, shapef, stridesf);
}


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

void atlas__FieldSet__add_fieldset(FieldSetImpl* This, FieldSetImpl* fieldset) {
    ATLAS_ASSERT(This != nullptr, "Reason: Use of uninitialised atlas_FieldSet");
    ATLAS_ASSERT(fieldset != nullptr, "Reason: Use of uninitialised atlas_FieldSet");
    for(int i = 0; i < fieldset->size(); i++) {
        This->add(fieldset->operator[](i));
    }
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
    fset.metadata() = metadata();
    return fset;
}

void FieldSet::set_dirty(bool value) const {
    get()->set_dirty(value);
}

//------------------------------------------------------------------------------------------------------

}  // namespace atlas
