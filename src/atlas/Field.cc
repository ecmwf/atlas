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

#include <eckit/exception/Exceptions.h>

#include "atlas/Field.h"
#include "atlas/FunctionSpace.h"

namespace atlas {

Field::Field(const std::string& name, const size_t nb_vars, FunctionSpace& function_space) :
  name_(name), nb_vars_(nb_vars), function_space_(function_space)
{
}

Field::~Field()
{
}

namespace {
template< typename DATA_TYPE >
DATA_TYPE* get_field_data( const Field& field )
{
  const FieldT<DATA_TYPE>* fieldT = dynamic_cast< const FieldT<DATA_TYPE>* >(&field);
  if( fieldT == NULL )
  {
    std::stringstream msg;
    msg << "Could not cast Field " << field.name()
        << " with data_type " << field.data_type() << " "
        << data_type_to_str<DATA_TYPE>();
    throw eckit::BadCast(msg.str(),Here());
  }
  return const_cast<DATA_TYPE*>(fieldT->data());
}
}

template <> const int*    Field::data<int   >() const { return get_field_data<int   >(*this); }
template <>       int*    Field::data<int   >()       { return get_field_data<int   >(*this); }
template <> const long*   Field::data<long  >() const { return get_field_data<long  >(*this); }
template <>       long*   Field::data<long  >()       { return get_field_data<long  >(*this); }
template <> const float*  Field::data<float >() const { return get_field_data<float >(*this); }
template <>       float*  Field::data<float >()       { return get_field_data<float >(*this); }
template <> const double* Field::data<double>() const { return get_field_data<double>(*this); }
template <>       double* Field::data<double>()       { return get_field_data<double>(*this); }


template<>
void FieldT<int>::halo_exchange() { function_space().halo_exchange(data_.data(),data_.size()); }
template<>
void FieldT<long>::halo_exchange() { function_space().halo_exchange(data_.data(),data_.size()); }
template<>
void FieldT<float>::halo_exchange() { function_space().halo_exchange(data_.data(),data_.size()); }
template<>
void FieldT<double>::halo_exchange() { function_space().halo_exchange(data_.data(),data_.size()); }

std::ostream& operator<<( std::ostream& os, const Field& f)
{
  f.print(os);
  return os;
}

// ------------------------------------------------------------------
// C wrapper interfaces to C++ routines

const char* atlas__Field__name (Field* This)
{
  return This->name().c_str();
}

const char* atlas__Field__data_type (Field* This)
{
  return This->data_type().c_str();
}

int atlas__Field__nb_vars (Field* This)
{
  return This->nb_vars();
}

Metadata* atlas__Field__metadata (Field* This)
{
  return &This->metadata();
}

FunctionSpace* atlas__Field__function_space (Field* This)
{
  return &This->function_space();
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

