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
#include "atlas/field/FieldT.h"
#include "atlas/field/FieldCreator.h"
#include "atlas/field/FieldTCreator.h"
#include "atlas/util/DataType.h"

#include "atlas/Grid.h"
#include "atlas/FunctionSpace.h"
#include "atlas/Mesh.h"

using atlas::field::FieldT;

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

  return NULL;
}

Field* Field::create(const std::string& name, const ArrayShape& shape, const std::string& datatype )
{
  eckit::ScopedPtr<field::FieldTCreator> creator
     (field::FieldTCreatorFactory::build("FieldT<"+datatype+">") );
  return creator->create_field(name,shape);
}

Field* Field::create( const ArrayShape& shape, const std::string& datatype )
{
  eckit::ScopedPtr<field::FieldTCreator> creator
     (field::FieldTCreatorFactory::build("FieldT<"+datatype+">") );
  return creator->create_field(std::string(),shape);
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

//Field::Field(DataType::kind_t kind, const ArrayShape& shape, const Indexing& index_type, const std::string& name):
//  name_(name)//, indexing_(index_type)
//{
//  array_.reset( ArrayBase::create(kind,shape) );
//}

Field::Field() :
  name_(), nb_levels_(0)
{
}

Field::Field(const std::string& name) :
  name_(name), nb_levels_(0)
{
}

Field::Field(const eckit::Parametrisation& params) :
  name_(), nb_levels_(0)
{
  params.get("name",name_);
}

Field::~Field()
{
}

//size_t Field::levels() const
//{
//  if( indexing_.level().size() == 1 )
//    return shape(indexing_.level());
//  else
//    return 0;
//}

//size_t Field::spectral() const
//{
//  if( indexing_.spectral().size() == 1 )
//    return shape(indexing_.spectral());
//  else
//    return 0;
//}

//namespace {
//template< typename DATA_TYPE >
//DATA_TYPE* get_field_data( const Field& field )
//{
//  try
//  {
//    return field.array_
//  }
//  catch ( eckit::BadCast& exception )
//  {

//  }

//  const FieldT<DATA_TYPE>* fieldT = dynamic_cast< const FieldT<DATA_TYPE>* >(&field);
//  if( fieldT == NULL )
//  {
//    std::stringstream msg;
//    msg << "Could not cast Field " << field.name()
//        << " with data_type " << field.datatype() << " to "
//        << DataType::datatype<DATA_TYPE>();
//    throw eckit::BadCast(msg.str(),Here());
//  }
//  return const_cast<DATA_TYPE*>(fieldT->data());
//}
//}

//template <> const int*    Field::data<int   >() const { return get_field_data<int   >(*this); }
//template <>       int*    Field::data<int   >()       { return get_field_data<int   >(*this); }
//template <> const long*   Field::data<long  >() const { return get_field_data<long  >(*this); }
//template <>       long*   Field::data<long  >()       { return get_field_data<long  >(*this); }
//template <> const float*  Field::data<float >() const { return get_field_data<float >(*this); }
//template <>       float*  Field::data<float >()       { return get_field_data<float >(*this); }
//template <> const double* Field::data<double>() const { return get_field_data<double>(*this); }
//template <>       double* Field::data<double>()       { return get_field_data<double>(*this); }

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

// ------------------------------------------------------------------
// C wrapper interfaces to C++ routines

Field* atlas__Field__create(eckit::Parametrisation* params)
{
  Field* field = 0;
  ATLAS_ERROR_HANDLING( field = Field::create(*params) );
  ASSERT(field);
  return field;
}

void atlas__Field__delete (Field* This)
{
  delete This;
}

const char* atlas__Field__name (Field* This)
{
  ATLAS_ERROR_HANDLING(
    ASSERT( This != NULL );
    return This->name().c_str();
  );
  return NULL;
}

void atlas__Field__datatype (Field* This, char* &datatype, int &size, int &allocated)
{
  ATLAS_ERROR_HANDLING(
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
    ASSERT( This != NULL );
    return This->size();
  );
  return 0;
}

int atlas__Field__rank (Field* This)
{
  ATLAS_ERROR_HANDLING(
    ASSERT( This != NULL );
    return This->rank();
  );
  return 0;
}

int atlas__Field__kind (Field* This)
{
  ATLAS_ERROR_HANDLING(
    ASSERT( This != NULL );
    return This->datatype().kind();
  );
  return 0;
}


double atlas__Field__bytes (Field* This)
{
  ATLAS_ERROR_HANDLING(
    ASSERT( This != NULL );
    return This->bytes();
  );
  return 0;
}

int atlas__Field__levels (Field* This)
{
  ATLAS_ERROR_HANDLING(
    ASSERT( This != NULL );
    return This->levels();
  );
  return 0;
}

Metadata* atlas__Field__metadata (Field* This)
{
  return &This->metadata();
}

int atlas__Field__has_functionspace(Field* This)
{
  ATLAS_ERROR_HANDLING(
    ASSERT( This != NULL );
    return This->has_functionspace();
  );
  return 0;
}

const char* atlas__Field__functionspace (Field* This)
{
  ATLAS_ERROR_HANDLING(
    ASSERT( This != NULL );
    return This->functionspace().c_str();
  );
  return NULL;
}


void atlas__Field__shapef (Field* This, int* &shape, int &rank)
{
  shape = const_cast<int*>(&This->shapef().front());
  rank = This->shapef().size();
}

void atlas__Field__data_shapef_int (Field* This, int* &field_data, int* &field_bounds, int &rank)
{
  field_data = &This->data<int>()[0];
  field_bounds = const_cast<int*>(&(This->shapef()[0]));
  rank = This->shapef().size();
}

void atlas__Field__data_shapef_long (Field* This, long* &field_data, int* &field_bounds, int &rank)
{
  field_data = &This->data<long>()[0];
  field_bounds = const_cast<int*>(&(This->shapef()[0]));
  rank = This->shapef().size();
}

void atlas__Field__data_shapef_float (Field* This, float* &field_data, int* &field_bounds, int &rank)
{
  field_data = &This->data<float>()[0];
  field_bounds = const_cast<int*>(&(This->shapef()[0]));
  rank = This->shapef().size();
}

void atlas__Field__data_shapef_double (Field* This, double* &field_data, int* &field_bounds, int &rank)
{
  field_data = &This->data<double>()[0];
  field_bounds = const_cast<int*>(&(This->shapef()[0]));
  rank = This->shapef().size();
}

// ------------------------------------------------------------------

} // namespace atlas

