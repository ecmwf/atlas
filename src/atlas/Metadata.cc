/*
 * (C) Copyright 1996-2015 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "atlas/Metadata.h"

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

namespace {
void throw_exception(const std::string& name)
{
  std::stringstream msg;
  msg << "Could not find metadata \"" << name << "\"";
  throw eckit::OutOfRange(msg.str(),Here());
}
}

Metadata& Metadata::set( const eckit::Properties& p )
{
  eckit::Properties::set(p);
  return *this;
}

template<> Metadata& Metadata::set(const std::string& name, const eckit::Properties& value)
{
  eckit::Properties::set(name,value);
  return *this;
}

template<> Metadata& Metadata::set(const std::string& name, const Metadata& value)
{
  eckit::Properties::set(name,value);
  return *this;
}

template<> Metadata& Metadata::set(const std::string& name, const std::vector<Metadata>& value)
{
  std::vector<eckit::Properties> properties_values(value.begin(),value.end());
  eckit::Properties::set(name, eckit::makeVectorValue(properties_values) );
  return *this;
}

template<> Metadata& Metadata::set(const std::string& name, const bool& value)
{
  eckit::Properties::set(name,value);
  return *this;
}

template<> Metadata& Metadata::set(const std::string& name, const int& value)
{
  eckit::Properties::set(name,value);
  return *this;
}

template<> Metadata& Metadata::set(const std::string& name, const long& value)
{
  eckit::Properties::set(name,value);
  return *this;
}

template<> Metadata& Metadata::set(const std::string& name, const double& value)
{
  eckit::Properties::set(name,value);
  return *this;
}

template<> Metadata& Metadata::set(const std::string& name, const size_t& value)
{
  eckit::Properties::set(name,value);
  return *this;
}

template<> Metadata& Metadata::set(const std::string& name, const std::string& value)
{
  eckit::Properties::set(name,value);
  return *this;
}

Metadata& Metadata::set(const std::string& name, const char* value)
{
  return set(name,std::string(value));
}


template<> Metadata& Metadata::set(const std::string& name, const std::vector<int>& value)
{
  eckit::Properties::set(name,eckit::makeVectorValue(value));
  return *this;
}

template<> Metadata& Metadata::set(const std::string& name, const std::vector<long>& value)
{
  eckit::Properties::set(name,eckit::makeVectorValue(value));
  return *this;
}

template<> Metadata& Metadata::set(const std::string& name, const std::vector<float>& value)
{
  eckit::Properties::set(name,eckit::makeVectorValue(value));
  return *this;
}

template<> Metadata& Metadata::set(const std::string& name, const std::vector<double>& value)
{
  eckit::Properties::set(name,eckit::makeVectorValue(value));
  return *this;
}

template<> Metadata& Metadata::set(const std::string& name, const std::vector<size_t>& value)
{
  eckit::Properties::set(name,eckit::makeVectorValue(value));
  return *this;
}

template<> Metadata& Metadata::set(const std::string& name, const Mesh& value)
{
  return set(name,value.id());
}

template<> Metadata& Metadata::set(const std::string& name, const Grid& value)
{
  return set(name,value.id());
}

template<> bool Metadata::get(const std::string& name) const
{
  if( !has(name) ) throw_exception(name);
  return eckit::Properties::get(name);
}

template<> int Metadata::get(const std::string& name) const
{
  if( !has(name) ) throw_exception(name);
  return eckit::Properties::get(name);
}

template<> long Metadata::get(const std::string& name) const
{
  if( !has(name) ) throw_exception(name);
  return eckit::Properties::get(name);
}

template<> float Metadata::get(const std::string& name) const
{
  if( !has(name) ) throw_exception(name);
  return double(eckit::Properties::get(name));
}

template<> double Metadata::get(const std::string& name) const
{
  if( !has(name) ) throw_exception(name);
  return eckit::Properties::get(name);
}

template<> size_t Metadata::get(const std::string& name) const
{
  if( !has(name) ) throw_exception(name);
  return eckit::Properties::get(name);
}

template<> std::string Metadata::get(const std::string& name) const
{
  if( !has(name) ) throw_exception(name);
  return eckit::Properties::get(name);
}

template<> std::vector<int> Metadata::get(const std::string& name) const
{
  if( !has(name) ) throw_exception(name);
  std::vector<eckit::Value> v = eckit::Properties::get(name);
  std::vector<int> value;
  value.assign(v.begin(),v.end());
  return value;
}

template<> std::vector<long> Metadata::get(const std::string& name) const
{
  if( !has(name) ) throw_exception(name);
  std::vector<eckit::Value> v = eckit::Properties::get(name);
  std::vector<long> value;
  value.assign(v.begin(),v.end());
  return value;
}

template<> std::vector<float> Metadata::get(const std::string& name) const
{
  if( !has(name) ) throw_exception(name);
  std::vector<eckit::Value> v = eckit::Properties::get(name);
  std::vector<float> value(v.size());
  for( size_t i=0; i<v.size(); ++i )
    value[i] = double( v[i] );
  return value;
}

template<> std::vector<double> Metadata::get(const std::string& name) const
{
  if( !has(name) ) throw_exception(name);
  std::vector<eckit::Value> v = eckit::Properties::get(name);
  std::vector<double> value;
  value.assign(v.begin(),v.end());
  return value;
}

template<> bool Metadata::get(const std::string& name, bool& value) const
{
  if( !has(name) ) return false;
  value = eckit::Properties::get(name);
  return true;
}

template<> bool Metadata::get(const std::string& name, int& value) const
{
  if( !has(name) ) return false;
  value = eckit::Properties::get(name);
  return true;
}

template<> bool Metadata::get(const std::string& name, long& value) const
{
  if( !has(name) ) return false;
  value = eckit::Properties::get(name);
  return true;
}

template<> bool Metadata::get(const std::string& name, float& value) const
{
  if( !has(name) ) return false;
  double v = eckit::Properties::get(name);
  value = v;
  return true;
}

template<> bool Metadata::get(const std::string& name, double& value) const
{
  if( !has(name) ) return false;
  value = eckit::Properties::get(name);
  return true;
}

template<> bool Metadata::get(const std::string& name, size_t& value) const
{
  if( !has(name) ) return false;
  value = eckit::Properties::get(name);
  return true;
}

template<> bool Metadata::get(const std::string& name, std::string& value) const
{
  if( !has(name) ) return false;
  value = std::string( eckit::Properties::get(name) );
  return true;
}

template<> bool Metadata::get(const std::string& name, std::vector<int>& value) const
{
  if(!has(name)) return false;
  std::vector<eckit::Value> v = operator[](name);
  value.assign(v.begin(),v.end());
  return true;
}
template<> bool Metadata::get(const std::string& name, std::vector<long>& value) const
{
  if(!has(name)) return false;
  std::vector<eckit::Value> v = operator[](name);
  value.assign(v.begin(),v.end());
  return true;
}
template<> bool Metadata::get(const std::string& name, std::vector<double>& value) const
{
  if(!has(name)) return false;
  std::vector<eckit::Value> v = operator[](name);
  value.assign(v.begin(),v.end());
  return true;
}

template<> bool Metadata::get(const std::string& name, std::vector<float>& value) const
{
  if(!has(name)) return false;
  std::vector<eckit::Value> v = operator[](name);
  value.resize(v.size());
  for( size_t i=0; i<v.size(); ++i )
    value[i] = double( v[i] );
  return true;
}

template<> Mesh& Metadata::get(const std::string& name) const
{
  if( !has(name) ) throw_exception(name);
  return Mesh::from_id(eckit::Properties::get(name));
}

template<> Grid& Metadata::get(const std::string& name) const
{
  if( !has(name) ) throw_exception(name);
  return Grid::from_id(eckit::Properties::get(name));
}

template<> eckit::Properties Metadata::get(const std::string& name) const
{
  if( !has(name) ) throw_exception(name);
  return eckit::Properties::get(name);
}

template<> std::vector<eckit::Properties> Metadata::get(const std::string& name) const
{
  if( !has(name) ) throw_exception(name);
  eckit::ValueList v = eckit::Properties::get(name);
  std::vector<eckit::Properties> values(v.size());
  for( size_t i=0; i<v.size(); ++i ) {
    values[i].set(v[i]);
  }
  return values;
}


bool Metadata::has(const std::string & name) const
{
  return Properties::has(name);
}


// ------------------------------------------------------------------
// C wrapper interfaces to C++ routines

Metadata* atlas__Metadata__new () {
  return new Metadata();
}

void atlas__Metadata__delete (Metadata* This) {
  delete This;
}

void atlas__Metadata__set_int (Metadata* This, const char* name, int value)
{
  ATLAS_ERROR_HANDLING( This->set( std::string(name), long(value) ) );
}
void atlas__Metadata__set_long (Metadata* This, const char* name, long value)
{
  ATLAS_ERROR_HANDLING( This->set( std::string(name), value ) );
}
void atlas__Metadata__set_float (Metadata* This, const char* name, float value)
{
  ATLAS_ERROR_HANDLING( This->set( std::string(name), double(value) ) );
}
void atlas__Metadata__set_double (Metadata* This, const char* name, double value)
{
  ATLAS_ERROR_HANDLING( This->set( std::string(name), value ) );
}
void atlas__Metadata__set_string (Metadata* This, const char* name, const char* value)
{
  ATLAS_ERROR_HANDLING( This->set( std::string(name), std::string(value) ) );
}
void atlas__Metadata__set_array_int (Metadata* This, const char* name, int value[], int size)
{
  ATLAS_ERROR_HANDLING(
    std::vector<int> v;
    v.assign(value,value+size);
    This->set( std::string(name), v );
  );
}
void atlas__Metadata__set_array_long (Metadata* This, const char* name, long value[], int size)
{
  ATLAS_ERROR_HANDLING(
    std::vector<long> v;
    v.assign(value,value+size);
    This->set( std::string(name), v );
  );
}
void atlas__Metadata__set_array_float (Metadata* This, const char* name, float value[], int size)
{
  ATLAS_ERROR_HANDLING(
    std::vector<float> v;
    v.assign(value,value+size);
    This->set( std::string(name), v );
  );
}
void atlas__Metadata__set_array_double (Metadata* This, const char* name, double value[], int size)
{
  ATLAS_ERROR_HANDLING(
    std::vector<double> v;
    v.assign(value,value+size);
    This->set( std::string(name), v );
  );
}
int atlas__Metadata__get_int (Metadata* This, const char* name)
{
  ATLAS_ERROR_HANDLING( return This->get<long>( std::string(name) ) );
  return 0;
}
long atlas__Metadata__get_long (Metadata* This, const char* name)
{
  ATLAS_ERROR_HANDLING( return This->get<long>( std::string(name) ) );
  return 0;
}
float atlas__Metadata__get_float (Metadata* This, const char* name)
{
  ATLAS_ERROR_HANDLING( return This->get<double>( std::string(name) ) );
  return 0;
}
double atlas__Metadata__get_double (Metadata* This, const char* name)
{
  ATLAS_ERROR_HANDLING( return This->get<double>( std::string(name) ) );
  return 0;
}
void atlas__Metadata__get_string( Metadata* This, const char* name, char* output_str, int max_len )
{
  ATLAS_ERROR_HANDLING(
    std::string s = This->get<std::string>( std::string(name) );
    if(s.size() > size_t(max_len))
    {
      std::stringstream msg;
      msg << "Cannot copy string `"<<s<<"` of metadata `"<<name<<"`"
             "in buffer of length " << max_len;
      throw eckit::OutOfRange(msg.str(),Here());
    }
    strcpy( output_str, s.c_str() );
    return
  );
  output_str = NULL;
}
void atlas__Metadata__get_array_int (Metadata* This, const char* name, int* &value, int& size, int& allocated)
{
  ATLAS_ERROR_HANDLING(
    std::vector<int> v = This->get< std::vector<int> >( std::string(name ) );
    size = v.size();
    value = new int[size];
    for( size_t j=0; j<v.size(); ++j ) value[j] = v[j];
    allocated = true;
  );
}
void atlas__Metadata__get_array_long (Metadata* This, const char* name, long* &value, int& size, int& allocated)
{
  ATLAS_ERROR_HANDLING(
    std::vector<long> v = This->get< std::vector<long> >( std::string(name ) );
    size = v.size();
    value = new long[size];
    for( size_t j=0; j<v.size(); ++j ) value[j] = v[j];
    allocated = true;
  );
}
void atlas__Metadata__get_array_float (Metadata* This, const char* name, float* &value, int& size, int& allocated)
{
  ATLAS_ERROR_HANDLING(
    std::vector<float> v = This->get< std::vector<float> >( std::string(name ) );
    size = v.size();
    value = new float[size];
    for( size_t j=0; j<v.size(); ++j ) value[j] = v[j];
    allocated = true;
  );
}
void atlas__Metadata__get_array_double (Metadata* This, const char* name, double* &value, int& size, int& allocated)
{
  ATLAS_ERROR_HANDLING(
    std::vector<double> v = This->get< std::vector<double> >( std::string(name ) );
    size = v.size();
    value = new double[size];
    for( size_t j=0; j<v.size(); ++j ) value[j] = v[j];
    allocated = true;
  );
}

void atlas__Metadata__set_grid (Metadata* This, const char* name, Grid* value)
{
  ATLAS_ERROR_HANDLING( This->set( std::string(name), *value ) );
}

void atlas__Metadata__set_mesh (Metadata* This, const char* name, Mesh* value)
{
  ATLAS_ERROR_HANDLING( This->set( std::string(name), *value ) );
}

Grid* atlas__Metadata__get_grid (Metadata* This, const char* name)
{
  Grid* value(NULL);
  ATLAS_ERROR_HANDLING( value = &This->get<Grid&>( std::string(name) ) );
  return value;
}

Mesh* atlas__Metadata__get_mesh (Metadata* This, const char* name)
{
  Mesh* value(NULL);
  ATLAS_ERROR_HANDLING( value = &This->get<Mesh&>( std::string(name) ) );
  return value;
}

int atlas__Metadata__has (Metadata* This, const char* name)
{
  ATLAS_ERROR_HANDLING( return This->has( std::string(name) ));
  return 0;
}

void atlas__Metadata__print (Metadata* This, std::ostream* channel)
{
  ATLAS_ERROR_HANDLING(
    ASSERT( This != NULL );
    ASSERT( channel != NULL );
    *channel << *This;
  );
}

void atlas__Metadata__json(Metadata* This, char* &json, int &size, int &allocated)
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
