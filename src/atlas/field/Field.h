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

#include <iosfwd>
#include <string>

#include "atlas/array/ArrayShape.h"
#include "atlas/array/DataType.h"
#include "atlas/array_fwd.h"
#include "atlas/util/ObjectHandle.h"

namespace eckit {
class Parametrisation;
}
namespace atlas {
namespace field {
class FieldImpl;
}
}  // namespace atlas
namespace atlas {
namespace util {
class Metadata;
}
}  // namespace atlas
namespace atlas {
class FunctionSpace;
}

namespace atlas {

class Field : public util::ObjectHandle<field::FieldImpl> {
public:
    using Handle::Handle;
    Field() = default;

    /// @brief Create field from parametrisation
    Field( const eckit::Parametrisation& );

    /// @brief Create field with given name, Datatype and ArrayShape
    Field( const std::string& name, array::DataType, const array::ArrayShape& = array::ArrayShape() );

    /// @brief Create field with given name, and take ownership of given Array
    Field( const std::string& name, array::Array* );

    /// @brief Create field by wrapping existing data, Datatype of template and
    /// ArraySpec
    template <typename DATATYPE>
    Field( const std::string& name, DATATYPE* data, const array::ArraySpec& );

    /// @brief Create field by wrapping existing data, Datatype of template and
    /// ArrayShape
    template <typename DATATYPE>
    Field( const std::string& name, DATATYPE* data, const array::ArrayShape& );

    // -- Conversion

    /// @brief Implicit conversion to Array
    operator const array::Array&() const;
    operator array::Array&();

    const array::Array& array() const;
    array::Array& array();

    bool valid() const { return get() != nullptr; }

    // -- Accessors

    /// @brief Access to raw data
    void* storage();

    /// @brief Internal data type of field
    array::DataType datatype() const;

    /// @brief Name associated to this field
    const std::string& name() const;

    /// @brief Rename this field
    void rename( const std::string& name );

    /// @brief Access to metadata associated to this field
    const util::Metadata& metadata() const;
    util::Metadata& metadata();

    /// @brief Resize field to given shape
    void resize( const array::ArrayShape& shape );

    void insert( idx_t idx1, idx_t size1 );

    /// @brief Shape of this field in Fortran style (reverse order of C style)
    const std::vector<int>& shapef() const;

    /// @brief Strides of this field in Fortran style (reverse order of C style)
    const std::vector<int>& stridesf() const;

    /// @brief Shape of this field (reverse order of Fortran style)
    const array::ArrayShape& shape() const;

    /// @brief Strides of this field
    const array::ArrayStrides& strides() const;

    /// @brief Shape of this field associated to index 'i'
    idx_t shape( idx_t i ) const;

    /// @brief Stride of this field associated to index 'i'
    idx_t stride( idx_t i ) const;

    /// @brief Number of values stored in this field
    idx_t size() const;

    /// @brief Rank of field
    idx_t rank() const;

    /// @brief Number of bytes occupied by the values of this field
    size_t bytes() const;

    bool contiguous() const;

    /// @brief Output information of field
    friend std::ostream& operator<<( std::ostream& os, const Field& v );

    /// @brief Output information of field plus raw data
    void dump( std::ostream& os ) const;

    /// Metadata that is more intrinsic to the Field, and queried often
    void set_levels( idx_t n );
    idx_t levels() const;

    /// Metadata that is more intrinsic to the Field, and queried often
    void set_variables( idx_t n );
    idx_t variables() const;

    void set_functionspace( const FunctionSpace& functionspace );
    const FunctionSpace& functionspace() const;

    /// @brief Return the memory footprint of the Field
    size_t footprint() const;

    bool dirty() const;

    void set_dirty( bool = true ) const;

    void haloExchange( bool on_device = false ) const;

    // -- dangerous methods
    template <typename DATATYPE>
    DATATYPE const* host_data() const;
    template <typename DATATYPE>
    DATATYPE* host_data();
    template <typename DATATYPE>
    DATATYPE const* device_data() const;
    template <typename DATATYPE>
    DATATYPE* device_data();
    template <typename DATATYPE>
    DATATYPE const* data() const;
    template <typename DATATYPE>
    DATATYPE* data();

    // -- Methods related to host-device synchronisation, requires
    // gridtools_storage
    void cloneToDevice() const;
    void cloneFromDevice() const;
    void syncHostDevice() const;
    bool hostNeedsUpdate() const;
    bool deviceNeedsUpdate() const;
    void reactivateDeviceWriteViews() const;
    void reactivateHostWriteViews() const;
};

extern template Field::Field( const std::string&, float*, const array::ArraySpec& );
extern template Field::Field( const std::string&, float*, const array::ArrayShape& );
extern template Field::Field( const std::string&, double*, const array::ArraySpec& );
extern template Field::Field( const std::string&, double*, const array::ArrayShape& );
extern template Field::Field( const std::string&, long*, const array::ArraySpec& );
extern template Field::Field( const std::string&, long*, const array::ArrayShape& );
extern template Field::Field( const std::string&, int*, const array::ArraySpec& );
extern template Field::Field( const std::string&, int*, const array::ArrayShape& );
extern template double const* Field::data() const;
extern template double* Field::data();
extern template float const* Field::data() const;
extern template float* Field::data();
extern template long const* Field::data() const;
extern template long* Field::data();
extern template int const* Field::data() const;
extern template int* Field::data();
extern template double const* Field::host_data() const;
extern template double* Field::host_data();
extern template float const* Field::host_data() const;
extern template float* Field::host_data();
extern template long const* Field::host_data() const;
extern template long* Field::host_data();
extern template int const* Field::host_data() const;
extern template int* Field::host_data();
extern template double const* Field::device_data() const;
extern template double* Field::device_data();
extern template float const* Field::device_data() const;
extern template float* Field::device_data();
extern template long const* Field::device_data() const;
extern template long* Field::device_data();
extern template int const* Field::device_data() const;
extern template int* Field::device_data();

//------------------------------------------------------------------------------------------------------

}  // namespace atlas
