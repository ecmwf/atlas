/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

/// @author Willem Deconinck
/// @date Sep 2014

#pragma once

#include <string>
#include <vector>

#include "eckit/memory/Owned.h"

#include "atlas/array.h"
#include "atlas/array/ArrayUtil.h"
#include "atlas/array/DataType.h"
#include "atlas/functionspace/FunctionSpace.h"
#include "atlas/util/Metadata.h"

namespace eckit {
class Parametrisation;
}

namespace atlas {
namespace field {

//----------------------------------------------------------------------------------------------------------------------

class FieldImpl : public eckit::Owned {
public:  // Static methods
    /// @brief Create field from parametrisation
    static FieldImpl* create( const eckit::Parametrisation& );

    /// @brief Create field with given name, Datatype and ArrayShape
    static FieldImpl* create( const std::string& name, array::DataType,
                              const array::ArrayShape& = array::ArrayShape() );

    /// @brief Create field with given name, Datatype of template and ArrayShape
    template <typename DATATYPE>
    static FieldImpl* create( const std::string& name, const array::ArrayShape& = array::ArrayShape() );

    /// @brief Create field with given name, and take ownership of given Array
    static FieldImpl* create( const std::string& name, array::Array* );

    /// @brief Create field with given name, and share ownership of given Array
    /// @note nawd: Not so sure we should go this route
    /// static FieldImpl* create( const std::string& name, const
    /// eckit::SharedPtr<Array>& );

    /// @brief Create field by wrapping existing data, Datatype of template and
    /// ArraySpec
    template <typename DATATYPE>
    static FieldImpl* wrap( const std::string& name, DATATYPE* data, const array::ArraySpec& );

    /// @brief Create field by wrapping existing data, Datatype of template and
    /// ArrayShape
    template <typename DATATYPE>
    static FieldImpl* wrap( const std::string& name, DATATYPE* data, const array::ArrayShape& );

private:  // Private constructors to force use of static create functions
    /// Allocate new Array internally
    FieldImpl( const std::string& name, array::DataType, const array::ArrayShape& );

    /// Transfer ownership of Array
    FieldImpl( const std::string& name, array::Array* );

    /// Share ownership of Array
    /// @note We could go this route...
    /// Field(const std::string& name, const eckit::SharedPtr<Array>& );

public:  // Destructor
    virtual ~FieldImpl();

    // -- Conversion

    /// @brief Implicit conversion to Array
    operator const array::Array&() const { return *array_; }
    operator array::Array&() { return *array_; }

    const array::Array& array() const { return *array_; }
    array::Array& array() { return *array_; }

    // -- Accessors

    /// @brief Access to raw data
    void* storage() { return array_->storage(); }

    /// @brief Internal data type of field
    array::DataType datatype() const { return array_->datatype(); }

    /// @brief Name associated to this field
    const std::string& name() const;

    /// @brief Rename this field
    void rename( const std::string& name ) { metadata().set( "name", name ); }

    /// @brief Access to metadata associated to this field
    const util::Metadata& metadata() const { return metadata_; }
    util::Metadata& metadata() { return metadata_; }

    /// @brief Resize field to given shape
    void resize( const array::ArrayShape& );

    void insert( size_t idx1, size_t size1 );

    /// @brief Shape of this field in Fortran style (reverse order of C style)
    const std::vector<int>& shapef() const { return array_->shapef(); }

    /// @brief Strides of this field in Fortran style (reverse order of C style)
    const std::vector<int>& stridesf() const { return array_->stridesf(); }

    /// @brief Shape of this field (reverse order of Fortran style)
    const array::ArrayShape& shape() const { return array_->shape(); }

    /// @brief Strides of this field
    const array::ArrayStrides& strides() const { return array_->strides(); }

    /// @brief Shape of this field associated to index 'i'
    size_t shape( size_t i ) const { return array_->shape( i ); }

    /// @brief Stride of this field associated to index 'i'
    size_t stride( size_t i ) const { return array_->stride( i ); }

    /// @brief Number of values stored in this field
    size_t size() const { return array_->size(); }

    /// @brief Rank of field
    size_t rank() const { return array_->rank(); }

    /// @brief Number of bytes occupied by the values of this field
    size_t bytes() const { return array_->bytes(); }

    /// @brief Output information of field
    friend std::ostream& operator<<( std::ostream& os, const FieldImpl& v );

    /// @brief Output information of field plus raw data
    void dump( std::ostream& os ) const;

    /// Metadata that is more intrinsic to the Field, and queried often
    void set_levels( size_t n ) { metadata().set( "levels", n ); }
    void set_variables( size_t n ) { metadata().set( "variables", n ); }
    size_t levels() const { return metadata().get<size_t>( "levels" ); }
    size_t variables() const { return metadata().get<size_t>( "variables" ); }

    void set_functionspace( const FunctionSpace& );
    const FunctionSpace& functionspace() const { return functionspace_; }

    /// @brief Return the memory footprint of the Field
    size_t footprint() const;

    // -- dangerous methods
    template <typename DATATYPE>
    DATATYPE const* host_data() const {
        return array_->host_data<DATATYPE>();
    }
    template <typename DATATYPE>
    DATATYPE* host_data() {
        return array_->host_data<DATATYPE>();
    }
    template <typename DATATYPE>
    DATATYPE const* device_data() const {
        return array_->device_data<DATATYPE>();
    }
    template <typename DATATYPE>
    DATATYPE* device_data() {
        return array_->device_data<DATATYPE>();
    }
    template <typename DATATYPE>
    DATATYPE const* data() const {
        return array_->host_data<DATATYPE>();
    }
    template <typename DATATYPE>
    DATATYPE* data() {
        return array_->host_data<DATATYPE>();
    }

    // -- Methods related to host-device synchronisation, requires
    // gridtools_storage
    void cloneToDevice() const { array_->cloneToDevice(); }
    void cloneFromDevice() const { array_->cloneFromDevice(); }
    void syncHostDevice() const { array_->syncHostDevice(); }
    bool hostNeedsUpdate() const { return array_->hostNeedsUpdate(); }
    bool deviceNeedsUpdate() const { return array_->deviceNeedsUpdate(); }
    void reactivateDeviceWriteViews() const { array_->reactivateDeviceWriteViews(); }
    void reactivateHostWriteViews() const { array_->reactivateHostWriteViews(); }

private:  // methods
    void print( std::ostream& os, bool dump = false ) const;

private:  // members
    mutable std::string name_;
    util::Metadata metadata_;
    array::Array* array_;
    FunctionSpace functionspace_;
};

//----------------------------------------------------------------------------------------------------------------------

template <typename DATATYPE>
FieldImpl* FieldImpl::create( const std::string& name, const array::ArrayShape& shape ) {
    return create( name, array::DataType::create<DATATYPE>(), shape );
}

template <typename DATATYPE>
FieldImpl* FieldImpl::wrap( const std::string& name, DATATYPE* data, const array::ArraySpec& spec ) {
    return create( name, array::Array::wrap( data, spec ) );
}

template <typename DATATYPE>
FieldImpl* FieldImpl::wrap( const std::string& name, DATATYPE* data, const array::ArrayShape& shape ) {
    return create( name, array::Array::wrap( data, shape ) );
}

//----------------------------------------------------------------------------------------------------------------------

// C wrapper interfaces to C++ routines
// #define Char char
extern "C" {
FieldImpl* atlas__Field__wrap_int_specf( const char* name, int data[], int rank, int shapef[], int stridesf[] );
FieldImpl* atlas__Field__wrap_long_specf( const char* name, long data[], int rank, int shapef[], int stridesf[] );
FieldImpl* atlas__Field__wrap_float_specf( const char* name, float data[], int rank, int shapef[], int stridesf[] );
FieldImpl* atlas__Field__wrap_double_specf( const char* name, double data[], int rank, int shapef[], int stridesf[] );
FieldImpl* atlas__Field__create( eckit::Parametrisation* params );
void atlas__Field__delete( FieldImpl* This );
const char* atlas__Field__name( FieldImpl* This );
void atlas__Field__datatype( FieldImpl* This, char*& datatype, int& size, int& allocated );
int atlas__Field__kind( FieldImpl* This );
int atlas__Field__rank( FieldImpl* This );
int atlas__Field__size( FieldImpl* This );
int atlas__Field__levels( FieldImpl* This );
double atlas__Field__bytes( FieldImpl* This );
void atlas__Field__shapef( FieldImpl* This, int*& shape, int& rank );
void atlas__Field__host_data_int_specf( FieldImpl* This, int*& field_data, int& rank, int*& field_shapef,
                                        int*& field_stridesf );
void atlas__Field__host_data_long_specf( FieldImpl* This, long*& field_data, int& rank, int*& field_shapef,
                                         int*& field_stridesf );
void atlas__Field__host_data_float_specf( FieldImpl* This, float*& field_data, int& rank, int*& field_shapef,
                                          int*& field_stridesf );
void atlas__Field__host_data_double_specf( FieldImpl* This, double*& field_data, int& rank, int*& field_shapef,
                                           int*& field_stridesf );
void atlas__Field__device_data_int_specf( FieldImpl* This, int*& field_data, int& rank, int*& field_shapef,
                                          int*& field_stridesf );
void atlas__Field__device_data_long_specf( FieldImpl* This, long*& field_data, int& rank, int*& field_shapef,
                                           int*& field_stridesf );
void atlas__Field__device_data_float_specf( FieldImpl* This, float*& field_data, int& rank, int*& field_shapef,
                                            int*& field_stridesf );
void atlas__Field__device_data_double_specf( FieldImpl* This, double*& field_data, int& rank, int*& field_shapef,
                                             int*& field_stridesf );
util::Metadata* atlas__Field__metadata( FieldImpl* This );
const functionspace::FunctionSpaceImpl* atlas__Field__functionspace( FieldImpl* This );
void atlas__Field__rename( FieldImpl* This, const char* name );
void atlas__Field__set_levels( FieldImpl* This, int levels );
void atlas__Field__set_functionspace( FieldImpl* This, const functionspace::FunctionSpaceImpl* functionspace );
int atlas__Field__host_needs_update( const FieldImpl* This );
int atlas__Field__device_needs_update( const FieldImpl* This );
void atlas__Field__clone_to_device( FieldImpl* This );
void atlas__Field__clone_from_device( FieldImpl* This );
void atlas__Field__sync_host_device( FieldImpl* This );
}
// #undef Char

//----------------------------------------------------------------------------------------------------------------------

}  // namespace field
}  // namespace atlas
