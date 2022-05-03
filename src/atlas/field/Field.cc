/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <iostream>
#include <string>

#include "atlas/library/config.h"

#include "atlas/field/Field.h"
#include "atlas/field/detail/FieldImpl.h"
#include "atlas/runtime/Exception.h"

#if ATLAS_HAVE_FUNCTIONSPACE
#include "atlas/functionspace/FunctionSpace.h"
#endif

namespace atlas {

// ------------------------------------------------------------------

std::ostream& operator<<(std::ostream& os, const Field& f) {
    os << (*f.get());
    return os;
}

Field::Field(const eckit::Parametrisation& config): Handle(Implementation::create(config)) {}

Field::Field(const std::string& name, array::DataType datatype, const array::ArrayShape& shape):
    Handle(Implementation::create(name, datatype, shape)) {}

Field::Field(const std::string& name, array::DataType datatype, array::ArraySpec&& spec):
    Handle(Implementation::create(name, datatype, std::move(spec))) {}

Field::Field(const std::string& name, array::Array* array): Handle(Implementation::create(name, array)) {}

template <typename DATATYPE>
Field::Field(const std::string& name, DATATYPE* data, const array::ArraySpec& spec):
    Handle(Implementation::wrap(name, data, spec)) {}

template <typename DATATYPE>
Field::Field(const std::string& name, DATATYPE* data, const array::ArrayShape& shape):
    Handle(Implementation::wrap(name, data, shape)) {}

/// @brief Implicit conversion to Array
Field::operator const array::Array&() const {
    return get()->array();
}
Field::operator array::Array&() {
    return get()->array();
}

const array::Array& Field::array() const {
    return get()->array();
}
array::Array& Field::array() {
    return get()->array();
}

/// @brief Clone
Field Field::clone(const eckit::Parametrisation& config) const {
    Field tmp(get()->name(), get()->datatype(), get()->shape());
    tmp.metadata() = this->metadata();
    tmp.set_functionspace(this->functionspace());
    array::Array::CopyPolicy cp;
      // To be set up via config. For now use default, as Array does not yet implement it.
    tmp.array().copy(this->array(),cp);
    return tmp;
}

// -- Accessors

/// @brief Access to raw data
void* Field::storage() {
    return get()->storage();
}

/// @brief Internal data type of field
array::DataType Field::datatype() const {
    return get()->datatype();
}

/// @brief Name associated to this field
const std::string& Field::name() const {
    return get()->name();
}

/// @brief Rename this field
void Field::rename(const std::string& name) {
    get()->rename(name);
}

/// @brief Access to metadata associated to this field
const util::Metadata& Field::metadata() const {
    return get()->metadata();
}
util::Metadata& Field::metadata() {
    return get()->metadata();
}

/// @brief Resize field to given shape
void Field::resize(const array::ArrayShape& shape) {
    get()->resize(shape);
}

void Field::insert(idx_t idx1, idx_t size1) {
    get()->insert(idx1, size1);
}

/// @brief Shape of this field in Fortran style (reverse order of C style)
const std::vector<int>& Field::shapef() const {
    return get()->shapef();
}

/// @brief Strides of this field in Fortran style (reverse order of C style)
const std::vector<int>& Field::stridesf() const {
    return get()->stridesf();
}

/// @brief Shape of this field (reverse order of Fortran style)
const array::ArrayShape& Field::shape() const {
    return get()->shape();
}

/// @brief Strides of this field
const array::ArrayStrides& Field::strides() const {
    return get()->strides();
}

/// @brief Shape of this field associated to index 'i'
idx_t Field::shape(idx_t i) const {
    return get()->shape(i);
}

/// @brief Stride of this field associated to index 'i'
idx_t Field::stride(idx_t i) const {
    return get()->stride(i);
}

/// @brief Number of values stored in this field
size_t Field::size() const {
    return get()->size();
}

/// @brief Rank of field
idx_t Field::rank() const {
    return get()->rank();
}

/// @brief Number of bytes occupied by the values of this field
size_t Field::bytes() const {
    return get()->bytes();
}

bool Field::contiguous() const {
    return array().contiguous();
}

/// @brief Output information of field plus raw data
void Field::dump(std::ostream& os) const {
    get()->dump(os);
}

/// Metadata that is more intrinsic to the Field, and queried often
void Field::set_levels(idx_t n) {
    get()->set_levels(n);
}
idx_t Field::levels() const {
    return get()->levels();
}

/// Metadata that is more intrinsic to the Field, and queried often
void Field::set_variables(idx_t n) {
    get()->set_variables(n);
}
idx_t Field::variables() const {
    return get()->variables();
}

void Field::set_horizontal_dimension(const std::vector<idx_t>& h_dim) {
    get()->set_horizontal_dimension(h_dim);
}

std::vector<idx_t> Field::horizontal_dimension() const {
    return get()->horizontal_dimension();
}

void Field::set_functionspace(const FunctionSpace& functionspace) {
#if ATLAS_HAVE_FUNCTIONSPACE
    get()->set_functionspace(functionspace);
#else
    throw_Exception("Atlas has been compiled without FunctionSpace support",Here());
#endif
}
const FunctionSpace& Field::functionspace() const {
#if ATLAS_HAVE_FUNCTIONSPACE
    return get()->functionspace();
#else
    throw_Exception("Atlas has been compiled without FunctionSpace support",Here());
#endif
}

/// @brief Return the memory footprint of the Field
size_t Field::footprint() const {
    return get()->footprint();
}

bool Field::dirty() const {
    return get()->dirty();
}

void Field::set_dirty(bool value) const {
    get()->set_dirty(value);
}

void Field::haloExchange(bool on_device) const {
    get()->haloExchange(on_device);
}

void Field::adjointHaloExchange(bool on_device) const {
    get()->adjointHaloExchange(on_device);
}

// -- Methods related to host-device synchronisation, requires gridtools_storage

void Field::updateHost() const {
    get()->updateHost();
}
void Field::updateDevice() const {
    get()->updateDevice();
}
void Field::syncHostDevice() const {
    get()->syncHostDevice();
}
bool Field::hostNeedsUpdate() const {
    return get()->hostNeedsUpdate();
}
bool Field::deviceNeedsUpdate() const {
    return get()->deviceNeedsUpdate();
}
void Field::setHostNeedsUpdate(bool v) const {
    return get()->setHostNeedsUpdate(v);
}
void Field::setDeviceNeedsUpdate(bool v) const {
    return get()->setDeviceNeedsUpdate(v);
}

// ------------------------------------------------------------------

template Field::Field(const std::string&, float*, const array::ArraySpec&);
template Field::Field(const std::string&, float*, const array::ArrayShape&);
template Field::Field(const std::string&, double*, const array::ArraySpec&);
template Field::Field(const std::string&, double*, const array::ArrayShape&);
template Field::Field(const std::string&, long*, const array::ArraySpec&);
template Field::Field(const std::string&, long*, const array::ArrayShape&);
template Field::Field(const std::string&, int*, const array::ArraySpec&);
template Field::Field(const std::string&, int*, const array::ArrayShape&);

// ------------------------------------------------------------------


}  // namespace atlas
