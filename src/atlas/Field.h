/*
 * (C) Copyright 1996-2015 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

/// @author Willem Deconinck
/// @author Tiago Quintino
/// @date Sep 2014

#ifndef atlas_Field_h
#define atlas_Field_h

#include <vector>
#include <string>

#include "eckit/memory/Owned.h"
#include "eckit/memory/SharedPtr.h"

#include "atlas/Config.h"
#include "atlas/Metadata.h"
#include "atlas/State.h"
#include "atlas/util/Array.h"

namespace eckit { class Parametrisation; }

namespace atlas {

//----------------------------------------------------------------------------------------------------------------------

class Field : public eckit::Owned {

public: // types

  typedef eckit::SharedPtr<Field> Ptr;

public: // Static methods

  static Field* create(const eckit::Parametrisation&);
  static Field* create(const ArrayShape&, const eckit::Parametrisation& = Config() );

public: // methods

// -- Constructor / Destructor
  Field();

  Field(const std::string& name);

  Field(const eckit::Parametrisation&);

  virtual ~Field();

// -- Accessors

  Ptr self() { return Ptr(this); }

  /// @brief Access to raw data
  template <typename DATA_TYPE> const DATA_TYPE* data() const  { return array_->data<DATA_TYPE>(); }
  template <typename DATA_TYPE>       DATA_TYPE* data()        { return array_->data<DATA_TYPE>(); }

  /// @brief Internal data type of field as string
  /// Any of [ int32 int64 real32 real64 ]
  std::string datatype() const { return array_->datatype(); }

  /// @brief Name associated to this field
  const std::string& name() const { return name_; }

  /// @brief Access to metadata associated to this field
  const Metadata& metadata() const { return metadata_; }
        Metadata& metadata()       { return metadata_; }

  void resize(const ArrayShape&);

  /// @brief Shape of this field in Fortran style (reverse order of C style)
  const std::vector<int>&    shapef()  const { return array_->shapef(); }

  /// @brief Shape of this field (reverse order of Fortran style)
  const ArrayShape& shape()   const { return array_->shape(); }

  /// @brief Strides of this field
  const ArrayStrides& strides() const { return array_->strides(); }

  /// @brief Shape of this field associated to index 'i'
  size_t shape (size_t i) const { return array_->shape(i); }

  /// @brief Stride of this field associated to index 'i'
  size_t stride(size_t i) const { return array_->stride(i); }

  /// @brief Number of values stored in this field
  size_t size() const { return array_->size(); }

  /// @brief Rank of field
  size_t rank() const { return array_->rank(); }

  /// @brief Number of bytes occupied by the values of this field
  virtual double bytes() const = 0;

  friend std::ostream& operator<<( std::ostream& os, const Field& v);

  virtual void dump(std::ostream& os) const = 0;

private: // methods

  virtual void print(std::ostream& os) const;

private: // friends

  // Allow a State::add() to change the name of the field, only if the field
  // has no name assigned, so it can be used for lookup later.
  friend Field& State::add( Field* );

private: // members

  std::string name_;

  Metadata metadata_;

protected: // members

  eckit::SharedPtr<ArrayBase> array_;

// End of class Field

//------------------------------------------------------------------------------------------------------

/***************************************************************************/
/* Public and private members and methods to be removed soon               */
/***************************************************************************/
private:

  const Grid* grid_;
  FunctionSpace* function_space_;

public:

  FunctionSpace& function_space() const;

  const Grid& grid() const;

  const Mesh& mesh() const;
        Mesh& mesh();

private:
  friend class FunctionSpace;
  void set_function_space(const FunctionSpace& function_space);

/***************************************************************************/
/* End Public and private members and  methods to be removed soon          */
/***************************************************************************/

};

//----------------------------------------------------------------------------------------------------------------------

#define Parametrisation eckit::Parametrisation
#define Char char
// C wrapper interfaces to C++ routines
extern "C"
{
  Field* atlas__Field__create(Parametrisation* params);
  void atlas__Field__delete (Field* This);
  const char* atlas__Field__name (Field* This);
  void atlas__Field__datatype (Field* This, Char* &datatype, int &size, int &allocated);
  int atlas__Field__rank (Field* This);
  int atlas__Field__size (Field* This);
  double atlas__Field__bytes (Field* This);
  void atlas__Field__shapef (Field* This, int* &shape, int &rank);
  void atlas__Field__data_shapef_int (Field* This, int* &field_data, int* &field_shapef, int &rank);
  void atlas__Field__data_shapef_long (Field* This, long* &field_data, int* &field_shapef, int &rank);
  void atlas__Field__data_shapef_float (Field* This, float* &field_data, int* &field_shapef, int &rank);
  void atlas__Field__data_shapef_double (Field* This, double* &field_data, int* &field_shapef, int &rank);
  Metadata* atlas__Field__metadata (Field* This);
  FunctionSpace* atlas__Field__function_space (Field* This);
}
#undef Parametrisation
#undef Char

//------------------------------------------------------------------------------------------------------

} // namespace atlas

#endif
