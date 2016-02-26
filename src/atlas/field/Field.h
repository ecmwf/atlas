/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

/// @author Willem Deconinck
/// @date Sep 2014

#ifndef atlas_Field_h
#define atlas_Field_h

#include <vector>
#include <string>
#include "eckit/memory/Owned.h"
#include "eckit/memory/SharedPtr.h"
#include "atlas/functionspace/FunctionSpace.h"
#include "atlas/util/DataType.h"
#include "atlas/util/array/ArrayUtil.h"
#include "atlas/util/array/Array.h"
#include "atlas/util/Metadata.h"

namespace eckit { class Parametrisation; }

namespace atlas {
namespace field {

//----------------------------------------------------------------------------------------------------------------------

class Field : public eckit::Owned {

public: // types

  typedef eckit::SharedPtr<Field> Ptr;

public: // Static methods

  /// @brief Create field from parametrisation
  static Field* create(const eckit::Parametrisation&);

  /// @brief Create field with given name, Datatype and ArrayShape
  static Field* create(
    const std::string& name, util::DataType,
    const util::array::ArrayShape& = util::array::ArrayShape());

  /// @brief Create field with given name, Datatype of template and ArrayShape
  template<typename DATATYPE>
  static Field* create( const std::string& name, const util::array::ArrayShape& = util::array::ArrayShape() );

  /// @brief Create field with given name, and take ownership of given Array
  static Field* create( const std::string& name, util::array::Array* );

  /// @brief Create field with given name, and share ownership of given Array
  /// @note nawd: Not so sure we should go this route
  /// static Field* create( const std::string& name, const eckit::SharedPtr<Array>& );

private: // Private constructors to force use of static create functions

  /// Allocate new Array internally
  Field(const std::string& name, util::DataType, const util::array::ArrayShape&);

  /// Transfer ownership of Array
  Field(const std::string& name, util::array::Array* );

  /// Share ownership of Array
  /// @note We could go this route...
  /// Field(const std::string& name, const eckit::SharedPtr<Array>& );

public: // Destructor
  virtual ~Field();

// -- Conversion

  /// @brief Implicit conversion to Array
  operator const util::array::Array&() const { return *array_; }
  operator util::array::Array&() { return *array_; }
  
  const util::array::Array& array() const { return *array_; }
        util::array::Array& array()       { return *array_; }

// -- Accessors

  /// @brief Access to raw data
  template <typename DATA_TYPE> const DATA_TYPE* data() const  { return array_->data<DATA_TYPE>(); }
  template <typename DATA_TYPE>       DATA_TYPE* data()        { return array_->data<DATA_TYPE>(); }

  /// @brief Internal data type of field
  util::DataType datatype() const { return array_->datatype(); }

  /// @brief Name associated to this field
  const std::string& name() const { return name_; }

  /// @brief Rename this field
  void rename(const std::string& name) { name_ = name; }

  /// @brief Access to metadata associated to this field
  const util::Metadata& metadata() const { return metadata_; }
        util::Metadata& metadata()       { return metadata_; }

  /// @brief Resize field to given shape
  void resize(const util::array::ArrayShape&);

  void insert(size_t idx1, size_t size1 );

  /// @brief Shape of this field in Fortran style (reverse order of C style)
  const std::vector<int>& shapef()  const { return array_->shapef(); }

  /// @brief Shape of this field (reverse order of Fortran style)
  const util::array::ArrayShape& shape() const { return array_->shape(); }

  /// @brief Strides of this field
  const util::array::ArrayStrides& strides() const { return array_->strides(); }

  /// @brief Shape of this field associated to index 'i'
  size_t shape (size_t i) const { return array_->shape(i); }

  /// @brief Stride of this field associated to index 'i'
  size_t stride(size_t i) const { return array_->stride(i); }

  /// @brief Number of values stored in this field
  size_t size() const { return array_->size(); }

  /// @brief Rank of field
  size_t rank() const { return array_->rank(); }

  /// @brief Number of bytes occupied by the values of this field
  double bytes() const { return array_->bytes(); }

  /// @brief Output information of field
  friend std::ostream& operator<<( std::ostream& os, const Field& v);

  /// @brief Output information of field plus raw data
  void dump(std::ostream& os) const;

  /// Metadata that is more intrinsic to the Field, and queried often
  bool has_levels() const { return nb_levels_!= 0; }
  void set_levels(size_t n) { nb_levels_ = n; }
  size_t levels() const { return std::max(1ul,nb_levels_); }

  void set_functionspace(const functionspace::FunctionSpace*);
  functionspace::FunctionSpace* functionspace() const { return functionspace_; }

private: // methods

  void print(std::ostream& os) const;

private: // members

  std::string name_;
  size_t nb_levels_;
  util::Metadata metadata_;
  util::array::Array* array_;
  functionspace::FunctionSpace* functionspace_;
};

//----------------------------------------------------------------------------------------------------------------------

template<typename DATATYPE>
Field* Field::create(
    const std::string& name,
    const util::array::ArrayShape& shape )
{
    return create(name, util::DataType::create<DATATYPE>(), shape);
}

//----------------------------------------------------------------------------------------------------------------------

// C wrapper interfaces to C++ routines
#define Parametrisation eckit::Parametrisation
#define functionspace_FunctionSpace functionspace::FunctionSpace
#define util_Metadata util::Metadata
#define Char char

extern "C"
{
  Field* atlas__Field__create(Parametrisation* params);
  void atlas__Field__delete (Field* This);
  const char* atlas__Field__name (Field* This);
  void atlas__Field__datatype (Field* This, Char* &datatype, int &size, int &allocated);
  int atlas__Field__kind (Field* This);
  int atlas__Field__rank (Field* This);
  int atlas__Field__size (Field* This);
  int atlas__Field__levels (Field* This);
  double atlas__Field__bytes (Field* This);
  void atlas__Field__shapef (Field* This, int* &shape, int &rank);
  void atlas__Field__data_shapef_int (Field* This, int* &field_data, int* &field_shapef, int &rank);
  void atlas__Field__data_shapef_long (Field* This, long* &field_data, int* &field_shapef, int &rank);
  void atlas__Field__data_shapef_float (Field* This, float* &field_data, int* &field_shapef, int &rank);
  void atlas__Field__data_shapef_double (Field* This, double* &field_data, int* &field_shapef, int &rank);
  util_Metadata* atlas__Field__metadata (Field* This);
  functionspace_FunctionSpace* atlas__Field__functionspace (Field* This);
}
#undef Parametrisation
#undef functionspace_FunctionSpace
#undef util_Metadata
#undef Char

//------------------------------------------------------------------------------------------------------

} // namespace field
} // namespace atlas

#endif
