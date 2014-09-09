/*
 * (C) Copyright 1996-2014 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */



#include <typeinfo>       // std::bad_cast
#include <sstream>
#include <stdexcept>

#include "atlas/mesh/Field.h"
#include "atlas/mesh/FunctionSpace.h"

namespace atlas {

Field::Field(const std::string& name, const int nb_vars, FunctionSpace& function_space) :
    name_(name), nb_vars_(nb_vars), function_space_(function_space)
{
}

Field::~Field()
{
}

template <>
int* Field::data<int>()
{
  try {
    return dynamic_cast< FieldT<int>& >(*this).data();
  } 
  catch (std::bad_cast& e) {
    std::stringstream msg;
    msg << "Could not cast Field " << name() 
        << " with data_type " << data_type() << " to int32";
    throw std::runtime_error(msg.str());
  }
}

template <>
int const* Field::data<int>() const
{
  try {
    return dynamic_cast< const FieldT<int>& >(*this).data();
  } 
  catch (std::bad_cast& e) {
    std::stringstream msg;
    msg << "Could not cast Field " << name() 
        << " with data_type " << data_type() << " to int32";
    throw std::runtime_error(msg.str());
  }
}

template <>
float* Field::data<float>()
{
  try {
    return dynamic_cast< FieldT<float>& >(*this).data();
  } 
  catch (std::bad_cast& e) {
    std::stringstream msg;
    msg << "Could not cast Field " << name() 
        << " with data_type " << data_type() << " to real32";
    throw std::runtime_error(msg.str());
  }
}

template <>
float const* Field::data<float>() const
{
  try {
    return dynamic_cast< const FieldT<float>& >(*this).data();
  } 
  catch (std::bad_cast& e) {
    std::stringstream msg;
    msg << "Could not cast Field " << name() 
        << " with data_type " << data_type() << " to real32";
    throw std::runtime_error(msg.str());
  }
}

template <>
double* Field::data<double>()
{
  try {
    return dynamic_cast< FieldT<double>& >(*this).data();
  } 
  catch (std::bad_cast& e) {
    std::stringstream msg;
    msg << "Could not cast Field " << name() 
        << " with data_type " << data_type() << " to real64";
    throw std::runtime_error(msg.str());
  }
}

template <>
double const* Field::data<double>() const
{
  try {
    return dynamic_cast< const FieldT<double>& >(*this).data();
  } 
  catch (std::bad_cast& e) {
    std::stringstream msg;
    msg << "Could not cast Field " << name() 
        << " with data_type " << data_type() << " to real64";
    throw std::runtime_error(msg.str());
  }
}


template<>
void FieldT<int>::halo_exchange() { function_space().halo_exchange(data_.data(),data_.size()); }
template<>
void FieldT<float>::halo_exchange() { function_space().halo_exchange(data_.data(),data_.size()); }
template<>
void FieldT<double>::halo_exchange() { function_space().halo_exchange(data_.data(),data_.size()); }

std::ostream& operator<<( std::ostream& os, const Field& f)
{
	f.print(os);
	return os;
}

#ifdef ECKIT_HAVE_GRIB

void Field::grib(Field::Grib *g)
{
	grib_.reset(g);
}

Field::Grib& Field::grib() const
{
	ASSERT( grib_ );
	return *grib_;
}

#endif

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


void atlas__Field__data_boundsf_double (Field* This, double* &field_data, int* &field_bounds, int &rank)
{
  field_data = &This->data<double>()[0];
  field_bounds = const_cast<int*>(&(This->boundsf()[0]));
  rank = This->boundsf().size();
}

void atlas__Field__data_boundsf_float (Field* This, float* &field_data, int* &field_bounds, int &rank)
{
  field_data = &This->data<float>()[0];
  field_bounds = const_cast<int*>(&(This->boundsf()[0]));
  rank = This->boundsf().size();
}

void atlas__Field__data_boundsf_int (Field* This, int* &field_data, int* &field_bounds, int &rank)
{
  field_data = &This->data<int>()[0];
  field_bounds = const_cast<int*>(&(This->boundsf()[0]));
  rank = This->boundsf().size();
}

// ------------------------------------------------------------------

} // namespace atlas

