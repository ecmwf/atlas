/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "atlas/Parametrisation.h"

#include <iostream>
#include <stdexcept>
#include <sstream>

#include "eckit/exception/Exceptions.h"
#include "eckit/parser/JSON.h"

#include "atlas/Mesh.h"
#include "atlas/Grid.h"
#include "atlas/FunctionSpace.h"
#include "atlas/runtime/ErrorHandling.h"

using std::string;

namespace atlas {

Parametrisation::Parametrisation(): inherited_(NULL) {}

Parametrisation::Parametrisation(const eckit::Properties& p): inherited_(NULL), delegate_(p) {}

Parametrisation::Parametrisation(const eckit::Parametrisation& p) : inherited_(NULL)
{
  if( const atlas::Parametrisation* params = dynamic_cast<const atlas::Parametrisation*>(&p) )
  {
    delegate_ = params->delegate_;
  }
  else
  {
    inherited_ = &p;
  }
}


Parametrisation& Parametrisation::set(const std::string& name, const char* value)
{
  delegate_.set(name,value);
  return *this;
}

//=================================================
// Functions overloading eckit::Parametrisation

bool Parametrisation::has(const std::string & name) const
{
  bool found(false);
  if( inherited_ )
    found = inherited_->has(name);
  found = found || delegate_.has(name);
  return found;
}

bool Parametrisation::get(const std::string& name, std::string& value) const
{
  bool found(false);
  if( inherited_ )
  {
    found = inherited_->get(name,value);
  }

  found = found || delegate_.get(name,value);

  return found;
}
bool Parametrisation::get(const std::string& name, bool& value) const
{
  bool found(false);
  if( inherited_ )
    found = inherited_->get(name,value);
  found = found || delegate_.get(name,value);
  return found;
}
bool Parametrisation::get(const std::string& name, int& value) const
{
  bool found(false);
  if( inherited_ )
    found = inherited_->get(name,value);
  found = found || delegate_.get(name,value);
  return found;
}
bool Parametrisation::get(const std::string& name, long& value) const
{
  bool found(false);
  if( inherited_ )
    found = inherited_->get(name,value);
  found = found || delegate_.get(name,value);
  return found;
}
bool Parametrisation::get(const std::string& name, size_t& value) const
{
  bool found(false);
  if( inherited_ )
    found = inherited_->get(name,value);
  found = found || delegate_.get(name,value);
  return found;
}
bool Parametrisation::get(const std::string& name, float& value) const
{
  bool found(false);
  if( inherited_ )
    found = inherited_->get(name,value);
  found = found || delegate_.get(name,value);
  return found;
}
bool Parametrisation::get(const std::string& name, double& value) const
{
  bool found(false);
  if( inherited_ )
    found = inherited_->get(name,value);
  found = found || delegate_.get(name,value);
  return found;
}
bool Parametrisation::get(const std::string& name, std::vector<int>& value) const
{
  bool found(false);
  if( inherited_ )
    found = inherited_->get(name,value);
  found = found || delegate_.get(name,value);
  return found;
}
bool Parametrisation::get(const std::string& name, std::vector<long>& value) const
{
  bool found(false);
  if( inherited_ )
    found = inherited_->get(name,value);
  found = found || delegate_.get(name,value);
  return found;
}
bool Parametrisation::get(const std::string& name, std::vector<float>& value) const
{
  bool found(false);
  if( inherited_ )
    found = inherited_->get(name,value);
  found = found || delegate_.get(name,value);
  return found;
}
bool Parametrisation::get(const std::string& name, std::vector<double>& value) const
{
  bool found(false);
  if( inherited_ )
    found = inherited_->get(name,value);
  found = found || delegate_.get(name,value);
  return found;
}

//==================================================================

// ------------------------------------------------------------------
// C wrapper interfaces to C++ routines

Parametrisation* atlas__Parametrisation__new () {
  return new Parametrisation();
}

void atlas__Parametrisation__delete (Parametrisation* This) {
  delete This;
}

void atlas__Parametrisation__set_int (Parametrisation* This, const char* name, int value)
{
  ATLAS_ERROR_HANDLING( This->set( std::string(name), long(value) ) );
}
void atlas__Parametrisation__set_long (Parametrisation* This, const char* name, long value)
{
  ATLAS_ERROR_HANDLING( This->set( std::string(name), value ) );
}
void atlas__Parametrisation__set_float (Parametrisation* This, const char* name, float value)
{
  ATLAS_ERROR_HANDLING( This->set( std::string(name), double(value) ) );
}
void atlas__Parametrisation__set_double (Parametrisation* This, const char* name, double value)
{
  ATLAS_ERROR_HANDLING( This->set( std::string(name), value ) );
}
void atlas__Parametrisation__set_string (Parametrisation* This, const char* name, const char* value)
{
  ATLAS_ERROR_HANDLING( This->set( std::string(name), std::string(value) ) );
}
void atlas__Parametrisation__set_array_int (Parametrisation* This, const char* name, int value[], int size)
{
  ATLAS_ERROR_HANDLING(
    std::vector<int> v;
    v.assign(value,value+size);
    This->set( std::string(name), v );
  );
}
void atlas__Parametrisation__set_array_long (Parametrisation* This, const char* name, long value[], int size)
{
  ATLAS_ERROR_HANDLING(
    std::vector<long> v;
    v.assign(value,value+size);
    This->set( std::string(name), v );
  );
}
void atlas__Parametrisation__set_array_float (Parametrisation* This, const char* name, float value[], int size)
{
  ATLAS_ERROR_HANDLING(
    std::vector<float> v;
    v.assign(value,value+size);
    This->set( std::string(name), v );
  );
}
void atlas__Parametrisation__set_array_double (Parametrisation* This, const char* name, double value[], int size)
{
  ATLAS_ERROR_HANDLING(
    std::vector<double> v;
    v.assign(value,value+size);
    This->set( std::string(name), v );
  );
}
int atlas__Parametrisation__get_int (Parametrisation* This, const char* name)
{
  int value;
  ATLAS_ERROR_HANDLING (
    if( ! This->get(std::string(name),value) )
        throw eckit::BadParameter("Parameter with name "+std::string(name)+" not found");
  );
  return value;
}
long atlas__Parametrisation__get_long (Parametrisation* This, const char* name)
{
  long value;
  ATLAS_ERROR_HANDLING (
    if( ! This->get(std::string(name),value) )
        throw eckit::BadParameter("Parameter with name "+std::string(name)+" not found");
  );
  return value;
}
float atlas__Parametrisation__get_float (Parametrisation* This, const char* name)
{
  float value;
  ATLAS_ERROR_HANDLING (
    if( ! This->get(std::string(name),value) )
        throw eckit::BadParameter("Parameter with name "+std::string(name)+" not found");
  );
  return value;
}
double atlas__Parametrisation__get_double (Parametrisation* This, const char* name)
{
  double value;
  ATLAS_ERROR_HANDLING (
    if( ! This->get(std::string(name),value) )
        throw eckit::BadParameter("Parameter with name "+std::string(name)+" not found");
  );
  return value;
}
void atlas__Parametrisation__get_string( Parametrisation* This, const char* name, char* output_str, int max_len )
{
  ATLAS_ERROR_HANDLING(
    std::string s;
    if( ! This->get(std::string(name),s) )
      throw eckit::BadParameter("Parameter with name "+std::string(name)+" not found");
    if( s.size() > max_len )
    {
      std::stringstream msg;
      msg << "Cannot copy string `"<<s<<"` of Parametrisation `"<<name<<"`"
             "in buffer of length " << max_len;
      throw eckit::OutOfRange(msg.str(),Here());
    }
    strcpy( output_str, s.c_str() );
    return
  );
  output_str = NULL;
}
void atlas__Parametrisation__get_array_int (Parametrisation* This, const char* name, int* &value, int& size, int& allocated)
{
  ATLAS_ERROR_HANDLING(
    std::vector<int> v;
    if( ! This->get(std::string(name),v) )
      throw eckit::BadParameter("Parameter with name "+std::string(name)+" not found");
    size = v.size();
    value = new int[size];
    for( size_t j=0; j<v.size(); ++j ) value[j] = v[j];
    allocated = true;
  );
}
void atlas__Parametrisation__get_array_long (Parametrisation* This, const char* name, long* &value, int& size, int& allocated)
{
  ATLAS_ERROR_HANDLING(
    std::vector<long> v;
    if( ! This->get(std::string(name),v) )
      throw eckit::BadParameter("Parameter with name "+std::string(name)+" not found");
    size = v.size();
    value = new long[size];
    for( size_t j=0; j<v.size(); ++j ) value[j] = v[j];
    allocated = true;
  );
}
void atlas__Parametrisation__get_array_float (Parametrisation* This, const char* name, float* &value, int& size, int& allocated)
{
  ATLAS_ERROR_HANDLING(
    std::vector<float> v;
    if( ! This->get(std::string(name),v) )
      throw eckit::BadParameter("Parameter with name "+std::string(name)+" not found");
    size = v.size();
    value = new float[size];
    for( size_t j=0; j<v.size(); ++j ) value[j] = v[j];
    allocated = true;
  );
}
void atlas__Parametrisation__get_array_double (Parametrisation* This, const char* name, double* &value, int& size, int& allocated)
{
  ATLAS_ERROR_HANDLING(
    std::vector<double> v;
    if( ! This->get(std::string(name),v) )
      throw eckit::BadParameter("Parameter with name "+std::string(name)+" not found");
    size = v.size();
    value = new double[size];
    for( size_t j=0; j<v.size(); ++j ) value[j] = v[j];
    allocated = true;
  );
}
int atlas__Parametrisation__has (Parametrisation* This, const char* name)
{
  ATLAS_ERROR_HANDLING( return This->has( std::string(name) ));
  return 0;
}

// ------------------------------------------------------------------

} // namespace atlas
