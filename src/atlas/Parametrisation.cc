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
#include "eckit/filesystem/PathName.h"
#include "eckit/parser/JSON.h"
#include "eckit/parser/JSONParser.h"

#include "atlas/Mesh.h"
#include "atlas/Grid.h"
#include "atlas/FunctionSpace.h"
#include "atlas/runtime/ErrorHandling.h"

using std::string;

namespace atlas {

Parametrisation::Parametrisation() {}

Parametrisation::Parametrisation(const eckit::Properties &p): delegate_(p) {}

Parametrisation::Parametrisation(std::istream& stream, const std::string &format )
{
  if( format != "json" ) {
    throw eckit::Exception("Not Implemented: Only json format is supported");
  }

  eckit::JSONParser parser( stream );
  set( eckit::Properties( parser.parse() ) );
}

Parametrisation::Parametrisation( const eckit::PathName& path )
{
  if( ! path.exists() ) {
    throw eckit::Exception("File "+std::string(path)+" does not exist.");
  }
  if( path.extension() == ".json" )
  {
    std::ifstream file(path.localPath());
    if (!file.is_open()) {
      throw eckit::Exception("Unable to open json file "+std::string(path),Here());
    }
    eckit::JSONParser parser( file );
    set( eckit::Properties( parser.parse() ) );
    file.close();
  }
  else
  {
    throw eckit::Exception("Only files with \".json\" extension supported for now. Found: "+path.extension());
  }
}



Parametrisation& Parametrisation::set(const eckit::Properties &p)
{
    delegate_.set(p);
    return *this;
}

Parametrisation& Parametrisation::set(const Parametrisation &p)
{
  delegate_.set(p.delegate_);
  return *this;
}


Parametrisation &Parametrisation::set(const std::string &name, const char *value) {
    delegate_.set(name, value);
    return *this;
}

Parametrisation& Parametrisation::set(const std::string& name, const Parametrisation& value )
{
  delegate_.set(name, value.delegate_);
  return *this;
}

Parametrisation& Parametrisation::set(const std::string& name, const std::vector<Parametrisation>& values )
{
  std::vector<Metadata> metadatavalues(values.size());
  for( size_t i=0; i<metadatavalues.size(); ++i )
    metadatavalues[i] = values[i].delegate_;
  delegate_.set(name, metadatavalues);
  return *this;
}



//=================================================
// Functions overloading eckit::Parametrisation

bool Parametrisation::has(const std::string &name) const {
    return delegate_.has(name);
}

template<class T>
bool Parametrisation::_get(const std::string &name, T &value) const {
    return delegate_.get(name, value);
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

bool Parametrisation::get(const std::string& name, Parametrisation& value) const {
  bool found = has(name);
  if( found ) {
    value.set( delegate_.get<eckit::Properties>(name) );
  }
  return found;
}
bool Parametrisation::get(const std::string& name, std::vector<Parametrisation>& value) const {
  bool found = has(name);
  if( found ) {
    std::vector<eckit::Properties> properties = delegate_.get< std::vector<eckit::Properties> >(name);
    value.resize(properties.size());
    for( size_t i=0; i<value.size(); ++i ) {
      value[i].set( properties[i] );
    }
  }
  return found;
}

eckit::JSON& operator<<(eckit::JSON& s, const Parametrisation& p)
{
  s << p.delegate_;
  return s;
}


//==================================================================

// ------------------------------------------------------------------
// C wrapper interfaces to C++ routines

Parametrisation* atlas__Parametrisation__new () {
  return new Parametrisation();
}

Parametrisation* atlas__Parametrisation__new_from_json (const char* json) {
  std::stringstream s;
  s << json;
  return new Parametrisation(s);
}


Parametrisation* atlas__Parametrisation__new_from_file (const char* path)
{
  return new Parametrisation( eckit::PathName(path) );
}

void atlas__Parametrisation__delete (Parametrisation* This) {
  ASSERT( This != 0 );
  delete This;
}

void atlas__Parametrisation__set_parametrisation (Parametrisation* This, const char* name, const Parametrisation* value)
{
  ATLAS_ERROR_HANDLING( This->set( std::string(name), *value ) );
}

void atlas__Parametrisation__set_parametrisation_list (Parametrisation* This, const char* name, const Parametrisation* value[], int size)
{
  std::vector<Parametrisation> params(size);
  for( size_t i=0; i<size; ++i )
  {
    params[i] = Parametrisation(*value[i]);
  }
  ATLAS_ERROR_HANDLING( This->set( std::string(name), params ) );
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

int atlas__Parametrisation__get_parametrisation (Parametrisation* This, const char* name, Parametrisation* value)
{
  ATLAS_ERROR_HANDLING ( if( ! This->get(std::string(name),*value) )  return false; );
  return true;
}

int atlas__Parametrisation__get_parametrisation_list (Parametrisation* This, const char* name, Parametrisation** &value, int &size, int &allocated)
{
  value = 0;
  ATLAS_ERROR_HANDLING (
    std::vector<Parametrisation> vector;
    if( ! This->get(std::string(name),vector) )  return false;
    size = vector.size();
    value = new Parametrisation*[size];
    allocated = true;
    for( size_t i=0; i<size; ++i ) {
      value[i] = new Parametrisation(vector[i]);
    }
  );
  return true;
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

void atlas__Parametrisation__json(Parametrisation* This, char* &json, int &size, int &allocated)
{
  std::stringstream s;
  eckit::JSON j(s);
  j.precision(16);
  j << *This;
  std::string json_str = s.str();
  size = json_str.size();
  json = new char[size+1]; allocated = true;
  strcpy(json,json_str.c_str());
  allocated = true;
}

// ------------------------------------------------------------------

} // namespace atlas
