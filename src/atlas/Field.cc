/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <typeinfo> // std::bad_cast
#include <sstream>
#include <stdexcept>

#include "eckit/exception/Exceptions.h"
#include "eckit/types/Types.h"

#include "atlas/runtime/ErrorHandling.h"
#include "atlas/FunctionSpace.h"
#include "atlas/Field.h"
#include "atlas/field/FieldCreator.h"
#include "atlas/util/DataType.h"

#include "atlas/Grid.h"
#include "atlas/FunctionSpace.h"
#include "atlas/Mesh.h"


namespace atlas {

template<class T>
inline std::ostream &operator<<(std::ostream &s, const std::vector<T> &v) {
    return eckit::__print_list(s, v);
}

// -------------------------------------------------------------------------
// Static functions

Field* Field::create(const eckit::Parametrisation& params)
{
  std::string creator_factory;
  if( params.get("creator",creator_factory) )
  {
    eckit::ScopedPtr<field::FieldCreator> creator
        (field::FieldCreatorFactory::build(creator_factory,params) );
    return creator->create_field(params);
  }
  else
    throw eckit::Exception("Could not find parameter 'creator' in Parametrisation for call to Field::create()");

  return 0;
}

Field* Field::create(DataType datatype, const ArrayShape& shape)
{
  return new Field(datatype,shape);
}

Field* Field::create(const std::string& name, DataType datatype, const ArrayShape& shape)
{
  return new Field(name,datatype,shape);
}

// -------------------------------------------------------------------------

Field::Field(DataType datatype, const ArrayShape& shape):
  name_(), nb_levels_(0)
{
  array_.reset( ArrayBase::create(datatype,shape) );
}

Field::Field(const std::string& name, DataType datatype, const ArrayShape& shape):
  name_(name), nb_levels_(0)
{
  array_.reset( ArrayBase::create(datatype,shape) );
}

Field::~Field()
{
}

void Field::dump(std::ostream& os) const
{
  print(os);
  array_->dump(os);
}

void Field::print(std::ostream& os) const
{
  os << "Field[name=" << name()
     << ",datatype=" << datatype().str()
     << ",size=" << size()
     << ",shape=" << shape()
     << ",strides=" << strides()
     << ",bytes=" << bytes()
     << ",metadata=" << metadata()
     << "]";
}

std::ostream& operator<<( std::ostream& os, const Field& f)
{
  f.print(os);
  return os;
}

void Field::resize(const ArrayShape& shape)
{
    array_->resize(shape);
}

template <> const float* Field::data() const
{
  try {
    return array_->data<float>();
  } catch( const eckit::BadCast& e ) {
    std::stringstream msg;
    msg << "Could not cast data for field " << name()
        << " with datatype " << array_->datatype().str() << " to "
        << DataType::str<float>();
    throw eckit::BadCast(msg.str(),Here());
  }
  return 0;
}
template <> float* Field::data()
{
  return const_cast<float*>(const_cast<const Field&>(*this).data<float>());
}

template <> const double* Field::data() const
{
  try {
    return array_->data<double>();
  } catch( const eckit::BadCast& e ) {
    std::stringstream msg;
    msg << "Could not cast data for field " << name()
        << " with datatype " << array_->datatype().str() << " to "
        << DataType::str<double>();
    throw eckit::BadCast(msg.str(),Here());
  }
  return 0;
}
template <> double* Field::data()
{
  return const_cast<double*>(const_cast<const Field&>(*this).data<double>());
}

template <> const int* Field::data() const
{
  try {
    return array_->data<int>();
  } catch( const eckit::BadCast& e ) {
    std::stringstream msg;
    msg << "Could not cast data for field " << name()
        << " with datatype " << array_->datatype().str() << " to "
        << DataType::str<int>();
    throw eckit::BadCast(msg.str(),Here());
  }
  return 0;
}
template <> int* Field::data()
{
  return const_cast<int*>(const_cast<const Field&>(*this).data<int>());
}

template <> const long* Field::data() const
{
  try {
    return array_->data<long>();
  } catch( const eckit::BadCast& e ) {
    std::stringstream msg;
    msg << "Could not cast data for field " << name()
        << " with datatype " << array_->datatype().str() << " to "
        << DataType::str<long>();
    throw eckit::BadCast(msg.str(),Here());
  }
  return 0;
}
template <> long* Field::data()
{
  return const_cast<long*>(const_cast<const Field&>(*this).data<long>());
}


// ------------------------------------------------------------------
// C wrapper interfaces to C++ routines

Field* atlas__Field__create(eckit::Parametrisation* params)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(params);
    Field* field = Field::create(*params);
    ASSERT(field);
    return field;
  );
  return 0;
}

void atlas__Field__delete (Field* This)
{
  delete This;
}

const char* atlas__Field__name (Field* This)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    return This->name().c_str();
  );
  return 0;
}

void atlas__Field__datatype (Field* This, char* &datatype, int &size, int &allocated)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    std::string s = This->datatype().str();
    size = s.size()+1;
    datatype = new char[size];
    strcpy(datatype,s.c_str());
    allocated = true;
  );
}

int atlas__Field__size (Field* This)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    return This->size();
  );
  return 0;
}

int atlas__Field__rank (Field* This)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    return This->rank();
  );
  return 0;
}

int atlas__Field__kind (Field* This)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    return This->datatype().kind();
  );
  return 0;
}


double atlas__Field__bytes (Field* This)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    return This->bytes();
  );
  return 0;
}

int atlas__Field__levels (Field* This)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    return This->levels();
  );
  return 0;
}

Metadata* atlas__Field__metadata (Field* This)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    return &This->metadata();
  );
  return 0;
}

int atlas__Field__has_functionspace(Field* This)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    return This->has_functionspace();
  );
  return 0;
}

const char* atlas__Field__functionspace (Field* This)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    return This->functionspace().c_str();
  );
  return 0;
}


void atlas__Field__shapef (Field* This, int* &shape, int &rank)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    shape = const_cast<int*>(&This->shapef().front());
    rank = This->shapef().size();
  );
}

void atlas__Field__data_shapef_int (Field* This, int* &field_data, int* &field_bounds, int &rank)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    field_data = &This->data<int>()[0];
    field_bounds = const_cast<int*>(&(This->shapef()[0]));
    rank = This->shapef().size();
  );
}

void atlas__Field__data_shapef_long (Field* This, long* &field_data, int* &field_bounds, int &rank)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    field_data = &This->data<long>()[0];
    field_bounds = const_cast<int*>(&(This->shapef()[0]));
    rank = This->shapef().size();
  );
}

void atlas__Field__data_shapef_float (Field* This, float* &field_data, int* &field_bounds, int &rank)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    field_data = &This->data<float>()[0];
    field_bounds = const_cast<int*>(&(This->shapef()[0]));
    rank = This->shapef().size();
  );
}

void atlas__Field__data_shapef_double (Field* This, double* &field_data, int* &field_bounds, int &rank)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    field_data = &This->data<double>()[0];
    field_bounds = const_cast<int*>(&(This->shapef()[0]));
    rank = This->shapef().size();
  );
}

// ------------------------------------------------------------------

} // namespace atlas

