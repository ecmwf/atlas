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

Parametrisation::Parametrisation(): inherited_(0) {}

Parametrisation::Parametrisation(const eckit::Properties &p): inherited_(0), delegate_(p) {}

Parametrisation::Parametrisation(const eckit::Parametrisation &p) : inherited_(0) {
    if ( const atlas::Parametrisation *params = dynamic_cast<const atlas::Parametrisation *>(&p) ) {
        delegate_ = params->delegate_;
    } else {
        inherited_ = &p;
    }
}


Parametrisation &Parametrisation::set(const std::string &name, const char *value) {
    delegate_.set(name, value);
    return *this;
}

//=================================================
// Functions overloading eckit::Parametrisation

bool Parametrisation::has(const std::string &name) const {
    bool found(false);
    if ( inherited_ )
        found = inherited_->has(name);
    found = found || delegate_.has(name);
    return found;
}

template<class T>
bool Parametrisation::_get(const std::string &name, T &value) const {
    bool found(false);
    if ( inherited_ ) {
        found = inherited_->get(name, value);
    }

    found = found || delegate_.get(name, value);

    return found;
}

bool Parametrisation::get(const std::string &name, std::string &value) const {
    return _get(name, value);
}

bool Parametrisation::get(const std::string &name, bool &value) const {
    return _get(name, value);
}

bool Parametrisation::get(const std::string &name, long &value) const {
    return _get(name, value);
}

bool Parametrisation::get(const std::string &name, size_t &value) const {
    return _get(name, value);
}

bool Parametrisation::get(const std::string &name, double &value) const {
    return _get(name, value);
}

bool Parametrisation::get(const std::string &name, std::vector<long> &value) const {
    return _get(name, value);
}

bool Parametrisation::get(const std::string &name, std::vector<double> &value) const {
    return _get(name, value);
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
int atlas__Parametrisation__get_int (Parametrisation* This, const char* name, int& value)
{
  long long_value;
  ATLAS_ERROR_HANDLING ( if( ! This->get(std::string(name),long_value) )  return false; );
  ASSERT( int(long_value) == long_value );
  value = long_value;
  return true;
}
int atlas__Parametrisation__get_long (Parametrisation* This, const char* name, long& value)
{
  ATLAS_ERROR_HANDLING ( if( ! This->get(std::string(name),value) )  return false; );
  return true;

}
int atlas__Parametrisation__get_float (Parametrisation* This, const char* name, float& value)
{
  double double_value;
  ATLAS_ERROR_HANDLING ( if ( ! This->get(std::string(name), double_value) ) return false ; );
  ASSERT(float(double_value) == double_value);
  value = double_value;
  return true;
}
int atlas__Parametrisation__get_double (Parametrisation* This, const char* name, double& value)
{
  ATLAS_ERROR_HANDLING ( if( ! This->get(std::string(name),value) )  return false; );
  return true;
}
int atlas__Parametrisation__get_string( Parametrisation* This, const char* name, char* &value, int &size, int &allocated )
{
  ATLAS_ERROR_HANDLING(
    std::string s;
    if( ! This->get(std::string(name),s) )
    {
      value = NULL;
      return false;
    }
    value = new char[s.size()+1];
    strcpy(value,s.c_str());
    allocated = true;
  );
  return true;
}
int atlas__Parametrisation__get_array_int (Parametrisation* This, const char* name, int* &value, int& size, int& allocated)
{
  ATLAS_ERROR_HANDLING(
    std::vector<long> v;
    if( ! This->get(std::string(name),v) )
      return false;
    size = v.size();
    value = new int[size];
    for ( size_t j = 0; j < v.size(); ++j ) {
      ASSERT(int(v[j]) == v[j]);
      value[j] = v[j];
    }
    allocated = true;
  );
  return true;
}
int atlas__Parametrisation__get_array_long (Parametrisation* This, const char* name, long* &value, int& size, int& allocated)
{
  ATLAS_ERROR_HANDLING(
    std::vector<long> v;
    if( ! This->get(std::string(name),v) )
      return false;
    size = v.size();
    value = new long[size];
    for( size_t j=0; j<v.size(); ++j ) value[j] = v[j];
    allocated = true;
  );
  return true;
}
int atlas__Parametrisation__get_array_float (Parametrisation* This, const char* name, float* &value, int& size, int& allocated)
{
  ATLAS_ERROR_HANDLING(
    std::vector<double> v;
    if( ! This->get(std::string(name),v) )
      return false;
    size = v.size();
    value = new float[size];
    for ( size_t j = 0; j < v.size(); ++j ) {
      ASSERT(float(v[j]) == v[j]);
      value[j] = v[j];
    }
    allocated = true;
  );
    return true;
}
int atlas__Parametrisation__get_array_double (Parametrisation* This, const char* name, double* &value, int& size, int& allocated)
{
  ATLAS_ERROR_HANDLING(
    std::vector<double> v;
    if( ! This->get(std::string(name),v) )
      return false;
    size = v.size();
    value = new double[size];
    for( size_t j=0; j<v.size(); ++j ) value[j] = v[j];
    allocated = true;
  );
  return true;
}

int atlas__Parametrisation__has (Parametrisation *This, const char *name) {
    ATLAS_ERROR_HANDLING( return This->has( std::string(name) ));
    return 0;
}

// ------------------------------------------------------------------

} // namespace atlas
