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

#include "atlas/field/Field.h"
#include "atlas/field/detail/FieldImpl.h"

namespace atlas {

// ------------------------------------------------------------------

std::ostream& operator<<( std::ostream& os, const Field& f ) {
    os << ( *f.field_ );
    return os;
}

Field::Field() : field_( nullptr ) {}

Field::Field( const Field& field ) : field_( field.field_ ) {
    field_->attach();
}

Field::Field( const Implementation* field ) : field_( const_cast<Implementation*>( field ) ) {
    field_->attach();
}

Field::Field( const eckit::Parametrisation& config ) : field_( Implementation::create( config ) ) {
    field_->attach();
}

Field::Field( const std::string& name, array::DataType datatype, const array::ArrayShape& shape ) :
    field_( Implementation::create( name, datatype, shape ) ) {
    field_->attach();
}

Field::Field( const std::string& name, array::Array* array ) : field_( Implementation::create( name, array ) ) {
    field_->attach();
}

template <>
Field::Field( const std::string& name, double* data, const array::ArraySpec& spec ) :
    field_( Implementation::wrap( name, data, spec ) ) {
    field_->attach();
}

template <>
Field::Field( const std::string& name, double* data, const array::ArrayShape& shape ) :
    field_( Implementation::wrap( name, data, shape ) ) {
    field_->attach();
}

template <>
Field::Field( const std::string& name, float* data, const array::ArraySpec& spec ) :
    field_( Implementation::wrap( name, data, spec ) ) {
    field_->attach();
}

template <>
Field::Field( const std::string& name, float* data, const array::ArrayShape& shape ) :
    field_( Implementation::wrap( name, data, shape ) ) {
    field_->attach();
}

template <>
Field::Field( const std::string& name, long* data, const array::ArraySpec& spec ) :
    field_( Implementation::wrap( name, data, spec ) ) {
    field_->attach();
}

template <>
Field::Field( const std::string& name, long* data, const array::ArrayShape& shape ) :
    field_( Implementation::wrap( name, data, shape ) ) {
    field_->attach();
}

template <>
Field::Field( const std::string& name, int* data, const array::ArraySpec& spec ) :
    field_( Implementation::wrap( name, data, spec ) ) {
    field_->attach();
}

template <>
Field::Field( const std::string& name, int* data, const array::ArrayShape& shape ) :
    field_( Implementation::wrap( name, data, shape ) ) {
    field_->attach();
}

Field::~Field() {
    if ( field_ ) {
        field_->detach();
        if ( not field_->owners() ) { delete field_; }
    }
}

const Field& Field::operator=( const Field& other ) {
    if ( field_ != other.field_ ) {
        if ( field_ ) {
            if ( not field_->owners() ) { delete field_; }
        }
        field_ = other.field_;
        field_->attach();
    }
    return *this;
}

/// @brief Implicit conversion to Array
Field::operator const array::Array&() const {
    return field_->array();
}
Field::operator array::Array&() {
    return field_->array();
}

const array::Array& Field::array() const {
    return field_->array();
}
array::Array& Field::array() {
    return field_->array();
}

// -- Accessors

/// @brief Access to raw data
void* Field::storage() {
    return field_->storage();
}

/// @brief Internal data type of field
array::DataType Field::datatype() const {
    return field_->datatype();
}

/// @brief Name associated to this field
const std::string& Field::name() const {
    return field_->name();
}

/// @brief Rename this field
void Field::rename( const std::string& name ) {
    field_->rename( name );
}

/// @brief Access to metadata associated to this field
const util::Metadata& Field::metadata() const {
    return field_->metadata();
}
util::Metadata& Field::metadata() {
    return field_->metadata();
}

/// @brief Resize field to given shape
void Field::resize( const array::ArrayShape& shape ) {
    field_->resize( shape );
}

void Field::insert( size_t idx1, size_t size1 ) {
    field_->insert( idx1, size1 );
}

/// @brief Shape of this field in Fortran style (reverse order of C style)
const std::vector<int>& Field::shapef() const {
    return field_->shapef();
}

/// @brief Strides of this field in Fortran style (reverse order of C style)
const std::vector<int>& Field::stridesf() const {
    return field_->stridesf();
}

/// @brief Shape of this field (reverse order of Fortran style)
const array::ArrayShape& Field::shape() const {
    return field_->shape();
}

/// @brief Strides of this field
const array::ArrayStrides& Field::strides() const {
    return field_->strides();
}

/// @brief Shape of this field associated to index 'i'
size_t Field::shape( size_t i ) const {
    return field_->shape( i );
}

/// @brief Stride of this field associated to index 'i'
size_t Field::stride( size_t i ) const {
    return field_->stride( i );
}

/// @brief Number of values stored in this field
size_t Field::size() const {
    return field_->size();
}

/// @brief Rank of field
size_t Field::rank() const {
    return field_->rank();
}

/// @brief Number of bytes occupied by the values of this field
size_t Field::bytes() const {
    return field_->bytes();
}

/// @brief Output information of field plus raw data
void Field::dump( std::ostream& os ) const {
    field_->dump( os );
}

/// Metadata that is more intrinsic to the Field, and queried often
void Field::set_levels( size_t n ) {
    field_->set_levels( n );
}
size_t Field::levels() const {
    return field_->levels();
}

/// Metadata that is more intrinsic to the Field, and queried often
void Field::set_variables( size_t n ) {
    field_->set_variables( n );
}
size_t Field::variables() const {
    return field_->variables();
}

void Field::set_functionspace( const FunctionSpace& functionspace ) {
    field_->set_functionspace( functionspace );
}
const FunctionSpace& Field::functionspace() const {
    return field_->functionspace();
}

/// @brief Return the memory footprint of the Field
size_t Field::footprint() const {
    return field_->footprint();
}

// -- dangerous methods
template <>
double const* Field::host_data() const {
    return field_->host_data<double>();
}
template <>
double* Field::host_data() {
    return field_->host_data<double>();
}
template <>
double const* Field::device_data() const {
    return field_->device_data<double>();
}
template <>
double* Field::device_data() {
    return field_->device_data<double>();
}
template <>
double const* Field::data() const {
    return field_->host_data<double>();
}
template <>
double* Field::data() {
    return field_->host_data<double>();
}

template <>
float const* Field::host_data() const {
    return field_->host_data<float>();
}
template <>
float* Field::host_data() {
    return field_->host_data<float>();
}
template <>
float const* Field::device_data() const {
    return field_->device_data<float>();
}
template <>
float* Field::device_data() {
    return field_->device_data<float>();
}
template <>
float const* Field::data() const {
    return field_->host_data<float>();
}
template <>
float* Field::data() {
    return field_->host_data<float>();
}

template <>
long const* Field::host_data() const {
    return field_->host_data<long>();
}
template <>
long* Field::host_data() {
    return field_->host_data<long>();
}
template <>
long const* Field::device_data() const {
    return field_->device_data<long>();
}
template <>
long* Field::device_data() {
    return field_->device_data<long>();
}
template <>
long const* Field::data() const {
    return field_->host_data<long>();
}
template <>
long* Field::data() {
    return field_->host_data<long>();
}

template <>
int const* Field::host_data() const {
    return field_->host_data<int>();
}
template <>
int* Field::host_data() {
    return field_->host_data<int>();
}
template <>
int const* Field::device_data() const {
    return field_->device_data<int>();
}
template <>
int* Field::device_data() {
    return field_->device_data<int>();
}
template <>
int const* Field::data() const {
    return field_->host_data<int>();
}
template <>
int* Field::data() {
    return field_->host_data<int>();
}

// -- Methods related to host-device synchronisation, requires gridtools_storage
void Field::cloneToDevice() const {
    field_->cloneToDevice();
}
void Field::cloneFromDevice() const {
    field_->cloneFromDevice();
}
void Field::syncHostDevice() const {
    field_->syncHostDevice();
}
bool Field::hostNeedsUpdate() const {
    return field_->hostNeedsUpdate();
}
bool Field::deviceNeedsUpdate() const {
    return field_->deviceNeedsUpdate();
}
void Field::reactivateDeviceWriteViews() const {
    field_->reactivateDeviceWriteViews();
}
void Field::reactivateHostWriteViews() const {
    field_->reactivateHostWriteViews();
}

// ------------------------------------------------------------------

}  // namespace atlas
