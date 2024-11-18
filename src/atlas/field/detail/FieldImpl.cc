/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <memory>
#include <sstream>

#include "atlas/library/config.h"

#include "atlas/array/MakeView.h"
#include "atlas/field/FieldCreator.h"
#include "atlas/field/detail/FieldImpl.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"

#if ATLAS_HAVE_FUNCTIONSPACE
#include "atlas/functionspace/FunctionSpace.h"
#endif

namespace atlas {
namespace field {

// -------------------------------------------------------------------------
// Static functions

FieldImpl* FieldImpl::create(const eckit::Parametrisation& params) {
    std::string creator_factory;
    if (params.get("creator", creator_factory)) {
        std::unique_ptr<field::FieldCreator> creator(field::FieldCreatorFactory::build(creator_factory, params));
        return creator->createField(params);
    }
    else {
        throw_Exception(
            "Could not find parameter 'creator' "
            "in Parametrisation for call to FieldImpl::create()");
    }
}

FieldImpl* FieldImpl::create(const std::string& name, array::DataType datatype, const array::ArrayShape& shape) {
    return new FieldImpl(name, datatype, shape);
}

FieldImpl* FieldImpl::create(const std::string& name, array::DataType datatype, array::ArraySpec&& spec) {
    return new FieldImpl(name, datatype, std::move(spec));
}

FieldImpl* FieldImpl::create(const std::string& name, array::Array* array) {
    return new FieldImpl(name, array);
}

// -------------------------------------------------------------------------

FieldImpl::FieldImpl(const std::string& name, array::DataType datatype, const array::ArrayShape& shape)
#if ATLAS_HAVE_FUNCTIONSPACE
    :functionspace_(new FunctionSpace())
#endif
{
    array_ = array::Array::create(datatype, shape);
    array_->attach();
    rename(name);
    set_levels(0);
    set_variables(0);
}

FieldImpl::FieldImpl(const std::string& name, array::DataType datatype, array::ArraySpec&& spec)
#if ATLAS_HAVE_FUNCTIONSPACE
    :functionspace_(new FunctionSpace())
#endif
{
    array_ = array::Array::create(datatype, std::move(spec));
    array_->attach();
    rename(name);
    set_levels(0);
    set_variables(0);
}


FieldImpl::FieldImpl(const std::string& name, array::Array* array)
#if ATLAS_HAVE_FUNCTIONSPACE
    :functionspace_(new FunctionSpace())
#endif
{
    array_ = array;
    array_->attach();
    rename(name);
    set_levels(0);
    set_variables(0);
}

FieldImpl::~FieldImpl() {
    for (FieldObserver* observer : field_observers_) {
        observer->onFieldDestruction(*this);
    }
    array_->detach();
    if (array_->owners() == 0) {
        for (auto& f : callback_on_destruction_) {
            f();
        }
        delete array_;
    }
#if ATLAS_HAVE_FUNCTIONSPACE
    delete functionspace_;
#endif
}

size_t FieldImpl::footprint() const {
    size_t size = sizeof(*this);
#if ATLAS_HAVE_FUNCTIONSPACE
    size += functionspace_->footprint();
#endif
    size += array_->footprint();
    size += metadata_.footprint();
    size += name_.capacity() * sizeof(std::string::value_type);
    return size;
}

bool FieldImpl::dirty() const {
    return metadata().getBool("dirty", true);
}

void FieldImpl::set_dirty(bool value) const {
    const_cast<FieldImpl&>(*this).metadata().set("dirty", value);
}

void FieldImpl::dump(std::ostream& os) const {
    print(os, true);
}

namespace {

template <typename T>
std::string vector_to_str(const std::vector<T>& t) {
    std::stringstream s;
    s << '[';
    for (size_t i = 0; i < t.size(); i++) {
        if (i != 0) {
            s << ',';
        }
        s << t[i];
    }
    s << ']';
    return s.str();
}

}  // namespace

void FieldImpl::rename(const std::string& name) {
    metadata().set("name", name);
    for (FieldObserver* observer : field_observers_) {
        observer->onFieldRename(*this);
    }
}

const std::string& FieldImpl::name() const {
    name_ = metadata().get<std::string>("name");
    return name_;
}

void FieldImpl::print(std::ostream& os, bool dump) const {
    os << "FieldImpl[name=" << name() << ",datatype=" << datatype().str() << ",size=" << size()
       << ",shape=" << vector_to_str(shape()) << ",strides=" << vector_to_str(strides())
#if !ATLAS_HAVE_GRIDTOOLS_STORAGE
       << ",bytes=" << bytes()
#endif
       << ",metadata=" << metadata();
    if (dump) {
        os << ",array=[";
        array_->dump(os);
        os << "]";
    }
    os << "]";
}

std::ostream& operator<<(std::ostream& os, const FieldImpl& f) {
    f.print(os);
    return os;
}

void FieldImpl::resize(const array::ArrayShape& shape) {
    array_->resize(shape);
}

void FieldImpl::insert(idx_t idx1, idx_t size1) {
    array_->insert(idx1, size1);
}

void FieldImpl::set_functionspace(const FunctionSpace& functionspace) {
#if ATLAS_HAVE_FUNCTIONSPACE
    *functionspace_ = functionspace;
#else
    throw_Exception("Atlas is compiled without FunctionSpace support", Here());
#endif
}

const FunctionSpace& FieldImpl::functionspace() const {
#if ATLAS_HAVE_FUNCTIONSPACE
    return *functionspace_;
#else
    throw_Exception("Atlas is compiled without FunctionSpace support", Here());
#endif
}

void FieldImpl::haloExchange(bool on_device) const {
    if (dirty()) {
#if ATLAS_HAVE_FUNCTIONSPACE
        ATLAS_ASSERT(functionspace());
        functionspace().haloExchange(Field(this), on_device);
        set_dirty(false);
#else
    throw_Exception("Atlas is compiled without FunctionSpace support", Here());
#endif
    }
}
void FieldImpl::adjointHaloExchange(bool on_device) const {
#if ATLAS_HAVE_FUNCTIONSPACE
    set_dirty();
    ATLAS_ASSERT(functionspace());
    functionspace().adjointHaloExchange(Field(this), on_device);
#else
    throw_Exception("Atlas is compiled without FunctionSpace support", Here());
#endif
}

void FieldImpl::attachObserver(FieldObserver& observer) const {
    if (std::find(field_observers_.begin(), field_observers_.end(), &observer) == field_observers_.end()) {
        field_observers_.push_back(&observer);
    }
}

void FieldImpl::detachObserver(FieldObserver& observer) const {
    field_observers_.erase(std::remove(field_observers_.begin(), field_observers_.end(), &observer),
                          field_observers_.end());
}


// ------------------------------------------------------------------

}  // namespace field
}  // namespace atlas
