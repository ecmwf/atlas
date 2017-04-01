/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <typeinfo>
#include <sstream>
#include <stdexcept>

#include "eckit/exception/Exceptions.h"
#include "eckit/memory/ScopedPtr.h"
#include "atlas/grid/Grid.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/field/Field.h"
#include "atlas/field/FieldCreator.h"
#include "atlas/functionspace/FunctionSpace.h"
#include "atlas/runtime/ErrorHandling.h"
#include "atlas/runtime/Log.h"
#include "atlas/array/MakeView.h"

namespace atlas {
namespace field {

// -------------------------------------------------------------------------
// Static functions

FieldImpl* FieldImpl::create(const eckit::Parametrisation& params)
{
  std::string creator_factory;
  if( params.get("creator",creator_factory) )
  {
    eckit::ScopedPtr<field::FieldCreator> creator
        (field::FieldCreatorFactory::build(creator_factory,params) );
    return creator->createField(params);
  }
  else
    throw eckit::Exception("Could not find parameter 'creator' "
                           "in Parametrisation for call to FieldImpl::create()");

  return 0;
}

FieldImpl* FieldImpl::create(
    const std::string& name,
    array::DataType           datatype,
    const array::ArrayShape&  shape)
{
  return new FieldImpl(name,datatype,shape);
}

FieldImpl* FieldImpl::create( const std::string& name, array::Array* array )
{
  return new FieldImpl(name,array);
}

// -------------------------------------------------------------------------

FieldImpl::FieldImpl(
    const std::string& name,
    array::DataType           datatype,
    const array::ArrayShape&  shape)
{
  array_ = array::Array::create(datatype,shape);
  array_->attach();
  rename(name);
  set_levels(0);
}


FieldImpl::FieldImpl(const std::string& name, array::Array* array)
{
  array_ = array;
  array_->attach();
  rename(name);
  set_levels(0);
}

FieldImpl::~FieldImpl()
{
  array_->detach();
  if( array_->owners() == 0 )
    delete array_;
}

size_t FieldImpl::footprint() const {
  size_t size = sizeof(*this);
  size += functionspace_.footprint();
  size += array_->footprint();
  size += metadata_.footprint();
  size += name_.capacity() * sizeof(std::string::value_type);
  return size;
}

void FieldImpl::dump(std::ostream& os) const
{
  print(os);
  array_->dump(os);
}

namespace {

template< typename T >
std::string vector_to_str(const std::vector<T>& t)
{
  std::stringstream s;
  s << '[';
  for(size_t i = 0; i < t.size(); i++) {
      if (i != 0)
          s << ',';
      s << t[i];
  }
  s << ']';
  return s.str();
}


}

const std::string& FieldImpl::name() const
{
  name_ = metadata().get<std::string>("name");
  return name_;
}

void FieldImpl::print(std::ostream& os) const
{
  os << "FieldImpl[name=" << name()
     << ",datatype=" << datatype().str()
     << ",size=" << size()
     << ",shape=" << vector_to_str( shape() )
     << ",strides=" << vector_to_str( strides() )
      #ifndef ATLAS_HAVE_GRIDTOOLS_STORAGE
     << ",bytes=" << bytes()
      #endif
     << ",metadata=" << metadata()
     << "]";
}

std::ostream& operator<<( std::ostream& os, const FieldImpl& f)
{
  f.print(os);
  return os;
}

void FieldImpl::resize(const array::ArrayShape& shape)
{
    array_->resize(shape);
}

void FieldImpl::insert(size_t idx1, size_t size1 )
{
    array_->insert(idx1,size1);
}


void FieldImpl::set_functionspace(const functionspace::FunctionSpace& functionspace)
{
  functionspace_ = functionspace;
}

// ------------------------------------------------------------------

Field::Field() :
  field_(nullptr) {
}

Field::Field( const Field& field ) :
  field_( field.field_ ) {
}

Field::Field( const FieldImpl* field ) : 
  field_( const_cast<FieldImpl*>(field) ) {
}

Field::Field(const eckit::Parametrisation& config) :
  field_( FieldImpl::create(config) ) {
} 

Field::Field(
  const std::string& name, array::DataType datatype,
  const array::ArrayShape& shape) :
  field_( FieldImpl::create(name, datatype, shape) ) {
}

Field::Field( const std::string& name, array::Array* array ) :
  field_( FieldImpl::create(name,array) ) {
}

// ------------------------------------------------------------------
// C wrapper interfaces to C++ routines

extern "C"
{

FieldImpl* atlas__Field__wrap_int_specf(const char* name, int data[], int rank, int shapef[], int stridesf[])
{
  ATLAS_ERROR_HANDLING(
    array::ArrayShape shape(rank);
    array::ArrayStrides strides(rank);
    size_t jf = rank-1;
    for( int j=0; j<rank; ++j )
    {
      shape  [j] = shapef  [jf];
      strides[j] = stridesf[jf];
      --jf;
    }
    FieldImpl* field;
    {
      Field wrapped(std::string(name),data,array::ArraySpec(shape,strides));
      field = wrapped.get();
      field->attach();
    }
    field->detach();
    ASSERT(field);
    return field;
  );
  return 0;
}

FieldImpl* atlas__Field__wrap_long_specf(const char* name, long data[], int rank, int shapef[], int stridesf[])
{
  ATLAS_ERROR_HANDLING(
    array::ArrayShape shape(rank);
    array::ArrayStrides strides(rank);
    size_t jf = rank-1;
    for( int j=0; j<rank; ++j )
    {
      shape  [j] = shapef  [jf];
      strides[j] = stridesf[jf];
      --jf;
    }
    FieldImpl* field;
    {
      Field wrapped(std::string(name),data,array::ArraySpec(shape,strides));
      field = wrapped.get();
      field->attach();
    }
    field->detach();
    ASSERT(field);
    return field;
  );
  return 0;
}

FieldImpl* atlas__Field__wrap_float_specf(const char* name, float data[], int rank, int shapef[], int stridesf[])
{
  ATLAS_ERROR_HANDLING(
    array::ArrayShape shape(rank);
    array::ArrayStrides strides(rank);
    size_t jf = rank-1;
    for( int j=0; j<rank; ++j )
    {
      shape  [j] = shapef  [jf];
      strides[j] = stridesf[jf];
      --jf;
    }
    FieldImpl* field;
    {
      Field wrapped(std::string(name),data,array::ArraySpec(shape,strides));
      field = wrapped.get();
      field->attach();
    }
    field->detach();
    ASSERT(field);
    return field;
  );
  return 0;
}

FieldImpl* atlas__Field__wrap_double_specf(const char* name, double data[], int rank, int shapef[], int stridesf[])
{
  ATLAS_ERROR_HANDLING(
    array::ArrayShape shape(rank);
    array::ArrayStrides strides(rank);
    size_t jf = rank-1;
    for( int j=0; j<rank; ++j )
    {
      shape  [j] = shapef  [jf];
      strides[j] = stridesf[jf];
      --jf;
    }
    FieldImpl* field;
    {
      Field wrapped(std::string(name),data,array::ArraySpec(shape,strides));
      field = wrapped.get();
      field->attach();
    }
    field->detach();
    ASSERT(field);
    return field;
  );
  return 0;
}

FieldImpl* atlas__Field__create(eckit::Parametrisation* params)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(params);
    FieldImpl* field;
    {
      Field f(*params);
      field = f.get();
      field->attach();
    }
    field->detach();
    
    ASSERT(field);
    return field;
  );
  return 0;
}

void atlas__Field__delete (FieldImpl* This)
{
  delete This;
}
  
const char* atlas__Field__name (FieldImpl* This)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    return This->name().c_str();
  );
  return 0;
}

void atlas__Field__datatype (FieldImpl* This, char* &datatype, int &size, int &allocated)
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

int atlas__Field__size (FieldImpl* This)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    return This->size();
  );
  return 0;
}

int atlas__Field__rank (FieldImpl* This)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    return This->rank();
  );
  return 0;
}

int atlas__Field__kind (FieldImpl* This)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    return This->datatype().kind();
  );
  return 0;
}


double atlas__Field__bytes (FieldImpl* This)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    return This->bytes();
  );
  return 0;
}

int atlas__Field__levels (FieldImpl* This)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    return This->levels();
  );
  return 0;
}

util::Metadata* atlas__Field__metadata (FieldImpl* This)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    return &This->metadata();
  );
  return 0;
}

int atlas__Field__has_functionspace(FieldImpl* This)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    return (This->functionspace() != 0);
  );
  return 0;
}

const functionspace::FunctionSpaceImpl* atlas__Field__functionspace (FieldImpl* This)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    return This->functionspace().get();
  );
  return 0;
}


void atlas__Field__shapef (FieldImpl* This, int* &shape, int &rank)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    shape = const_cast<int*>(&This->shapef().front());
    rank = This->shapef().size();
  );
}

void atlas__Field__host_data_int_specf (FieldImpl* This, int* &data, int &rank, int* &shapef, int* &stridesf)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    data = This->host_data<int>();
    shapef = const_cast<int*>(This->shapef().data());
    stridesf = const_cast<int*>(This->stridesf().data());
    rank = This->shapef().size();
  );
}

void atlas__Field__host_data_long_specf (FieldImpl* This, long* &data, int &rank, int* &shapef, int* &stridesf)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    data = This->host_data<long>();
    shapef = const_cast<int*>(This->shapef().data());
    stridesf = const_cast<int*>(This->stridesf().data());
    rank = This->shapef().size();
  );
}

void atlas__Field__host_data_float_specf (FieldImpl* This, float* &data, int &rank, int* &shapef, int* &stridesf)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    data = This->host_data<float>();
    shapef = const_cast<int*>(This->shapef().data());
    stridesf = const_cast<int*>(This->stridesf().data());
    rank = This->shapef().size();
  );
}

void atlas__Field__host_data_double_specf (FieldImpl* This, double* &data, int &rank, int* &shapef, int* &stridesf)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    data = This->host_data<double>();
    shapef = const_cast<int*>(This->shapef().data());
    stridesf = const_cast<int*>(This->stridesf().data());
    rank = This->shapef().size();
  );
}

void atlas__Field__device_data_int_specf (FieldImpl* This, int* &data, int &rank, int* &shapef, int* &stridesf)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    data = This->device_data<int>();
    shapef = const_cast<int*>(This->shapef().data());
    stridesf = const_cast<int*>(This->stridesf().data());
    rank = This->shapef().size();
  );
}

void atlas__Field__device_data_long_specf (FieldImpl* This, long* &data, int &rank, int* &shapef, int* &stridesf)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    data = This->device_data<long>();
    shapef = const_cast<int*>(This->shapef().data());
    stridesf = const_cast<int*>(This->stridesf().data());
    rank = This->shapef().size();
  );
}

void atlas__Field__device_data_float_specf (FieldImpl* This, float* &data, int &rank, int* &shapef, int* &stridesf)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    data = This->device_data<float>();
    shapef = const_cast<int*>(This->shapef().data());
    stridesf = const_cast<int*>(This->stridesf().data());
    rank = This->shapef().size();
  );
}

void atlas__Field__device_data_double_specf (FieldImpl* This, double* &data, int &rank, int* &shapef, int* &stridesf)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    data = This->device_data<double>();
    shapef = const_cast<int*>(This->shapef().data());
    stridesf = const_cast<int*>(This->stridesf().data());
    rank = This->shapef().size();
  );
}

int atlas__Field__is_on_host(const FieldImpl* This)
{
  return This->isOnHost();
}

int atlas__Field__is_on_device(const FieldImpl* This)
{
  return This->isOnDevice();
}

void atlas__Field__rename(FieldImpl* This, const char* name)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    This->rename( std::string(name) );
  );
}

void atlas__Field__set_levels(FieldImpl* This, int levels)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    This->set_levels(levels);
  );
}

void atlas__Field__set_functionspace(FieldImpl* This, const functionspace::FunctionSpaceImpl* functionspace)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    ASSERT(functionspace);
    This->set_functionspace( functionspace );
  );
}

void atlas__Field__clone_to_device(FieldImpl* This)
{
  This->cloneToDevice();
}
void atlas__Field__clone_from_device(FieldImpl* This)
{
  This->cloneFromDevice();
}


}
// ------------------------------------------------------------------

} // namespace field
} // namespace atlas

