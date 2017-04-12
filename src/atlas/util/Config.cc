/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <iostream>
#include <stdexcept>
#include <sstream>
#include <cstdarg>
#include "eckit/exception/Exceptions.h"
#include "eckit/filesystem/PathName.h"
#include "eckit/parser/JSON.h"
#include "eckit/parser/JSONParser.h"
#include "atlas/grid/Grid.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/util/Config.h"
#include "atlas/runtime/ErrorHandling.h"

using std::string;

namespace atlas {
namespace util {

Config::Config() {}

Config::Config(const eckit::Properties &p): delegate_(p) {}

Config::Config(std::istream& stream, const std::string &format )
{
  if( format != "json" ) {
    throw eckit::Exception("Not Implemented: Only json format is supported");
  }

  eckit::JSONParser parser( stream );
  set( eckit::Properties( parser.parse() ) );
}

Config::Config( const eckit::PathName& path )
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

// Config Config::operator&&(const Config& other) const
// {
//    Config config;
//    config.set(*this);
//    config.set(other);
//    return config;
// }

Config Config::operator|(const Config& other) const
{
   Config config;
   config.set(*this);
   config.set(other);
   return config;
}

Config& Config::set(const eckit::Properties &p)
{
    delegate_.set(p);
    return *this;
}

Config& Config::set(const Config &p)
{
  delegate_.set(p.delegate_);
  return *this;
}


Config &Config::set(const std::string &name, const char *value) {
    delegate_.set(name, value);
    return *this;
}

Config& Config::set(const std::string& name, const Config& value )
{
  delegate_.set(name, value.delegate_);
  return *this;
}

Config& Config::set(const std::string& name, const std::vector<Config>& values )
{
  std::vector<Metadata> metadatavalues(values.size());
  for( size_t i=0; i<metadatavalues.size(); ++i )
    metadatavalues[i] = values[i].delegate_;
  delegate_.set(name, metadatavalues);
  return *this;
}



//=================================================
// Functions overloading eckit::Parametrisation

bool Config::has(const std::string &name) const {
    return delegate_.has(name);
}

template<class T>
bool Config::_get(const std::string &name, T &value) const {
    return delegate_.get(name, value);
}

bool Config::get(const std::string &name, std::string &value) const {
    return _get(name, value);
}

bool Config::get(const std::string &name, bool &value) const {
    return _get(name, value);
}

bool Config::get(const std::string &name, long &value) const {
    return _get(name, value);
}

bool Config::get(const std::string &name, size_t &value) const {
    return _get(name, value);
}

bool Config::get(const std::string &name, double &value) const {
    return _get(name, value);
}

bool Config::get(const std::string &name, std::vector<long> &value) const {
    return _get(name, value);
}

bool Config::get(const std::string &name, std::vector<double> &value) const {
    return _get(name, value);
}

bool Config::get(const std::string& name, Config& value) const {
  bool found = has(name);
  if( found ) {
    value.set( delegate_.get<eckit::Properties>(name) );
  }
  return found;
}
bool Config::get(const std::string& name, std::vector<Config>& value) const {
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

eckit::JSON& operator<<(eckit::JSON& s, const Config& p)
{
  s << p.delegate_;
  return s;
}

std::ostream& operator<<(std::ostream& s, const Config& p)
{
  s << p.delegate_;
  return s;
}



//==================================================================

// ------------------------------------------------------------------
// C wrapper interfaces to C++ routines

Config* atlas__Config__new () {
  return new Config();
}

Config* atlas__Config__new_from_json (const char* json) {
  std::stringstream s;
  s << json;
  return new Config(s);
}


Config* atlas__Config__new_from_file (const char* path)
{
  return new Config( eckit::PathName(path) );
}

void atlas__Config__delete (Config* This) {
  ASSERT( This != 0 );
  delete This;
}

void atlas__Config__set_config (Config* This, const char* name, const Config* value)
{
  ATLAS_ERROR_HANDLING( This->set( std::string(name), *value ) );
}

void atlas__Config__set_config_list (Config* This, const char* name, const Config* value[], int size)
{
  std::vector<Config> params(size);
  for(int i = 0; i < size; ++i)
  {
    params[i] = Config(*value[i]);
  }
  ATLAS_ERROR_HANDLING( This->set( std::string(name), params ) );
}

void atlas__Config__set_int (Config* This, const char* name, int value)
{
  ATLAS_ERROR_HANDLING( This->set( std::string(name), long(value) ) );
}
void atlas__Config__set_long (Config* This, const char* name, long value)
{
  ATLAS_ERROR_HANDLING( This->set( std::string(name), value ) );
}
void atlas__Config__set_float (Config* This, const char* name, float value)
{
  ATLAS_ERROR_HANDLING( This->set( std::string(name), double(value) ) );
}
void atlas__Config__set_double (Config* This, const char* name, double value)
{
  ATLAS_ERROR_HANDLING( This->set( std::string(name), value ) );
}
void atlas__Config__set_string (Config* This, const char* name, const char* value)
{
  ATLAS_ERROR_HANDLING( This->set( std::string(name), std::string(value) ) );
}
void atlas__Config__set_array_int (Config* This, const char* name, int value[], int size)
{
  ATLAS_ERROR_HANDLING(
    std::vector<int> v;
    v.assign(value,value+size);
    This->set( std::string(name), v );
  );
}
void atlas__Config__set_array_long (Config* This, const char* name, long value[], int size)
{
  ATLAS_ERROR_HANDLING(
    std::vector<long> v;
    v.assign(value,value+size);
    This->set( std::string(name), v );
  );
}
void atlas__Config__set_array_float (Config* This, const char* name, float value[], int size)
{
  ATLAS_ERROR_HANDLING(
    std::vector<float> v;
    v.assign(value,value+size);
    This->set( std::string(name), v );
  );
}
void atlas__Config__set_array_double (Config* This, const char* name, double value[], int size)
{
  ATLAS_ERROR_HANDLING(
    std::vector<double> v;
    v.assign(value,value+size);
    This->set( std::string(name), v );
  );
}

int atlas__Config__get_config (Config* This, const char* name, Config* value)
{
  ATLAS_ERROR_HANDLING ( if( ! This->get(std::string(name),*value) )  return false; );
  return true;
}

int atlas__Config__get_config_list (Config* This, const char* name, Config** &value, int &size, int &allocated)
{
  value = 0;
  ATLAS_ERROR_HANDLING (
    std::vector<Config> vector;
    if( ! This->get(std::string(name),vector) )  return false;
    size = vector.size();
    value = new Config*[size];
    allocated = true;
    for(int i = 0; i < size; ++i) {
      value[i] = new Config(vector[i]);
    }
  );
  return true;
}

int atlas__Config__get_int (Config* This, const char* name, int& value)
{
  long long_value = value;
  ATLAS_ERROR_HANDLING ( if( ! This->get(std::string(name),long_value) )  return false; );
  ASSERT( int(long_value) == long_value );
  value = long_value;
  return true;
}
int atlas__Config__get_long (Config* This, const char* name, long& value)
{
  ATLAS_ERROR_HANDLING ( if( ! This->get(std::string(name),value) )  return false; );
  return true;

}
int atlas__Config__get_float (Config* This, const char* name, float& value)
{
  double double_value;
  ATLAS_ERROR_HANDLING ( if ( ! This->get(std::string(name), double_value) ) return false ; );
  value = double_value;
  return true;
}
int atlas__Config__get_double (Config* This, const char* name, double& value)
{
  ATLAS_ERROR_HANDLING ( if( ! This->get(std::string(name),value) )  return false; );
  return true;
}
int atlas__Config__get_string( Config* This, const char* name, char* &value, int &size, int &allocated )
{
  ATLAS_ERROR_HANDLING(
    std::string s;
    if( ! This->get(std::string(name),s) )
    {
      value = NULL;
      return false;
    }
    size = s.size()+1;
    value = new char[size];
    strcpy(value,s.c_str());
    allocated = true;
  );
  return true;
}
int atlas__Config__get_array_int (Config* This, const char* name, int* &value, int& size, int& allocated)
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
int atlas__Config__get_array_long (Config* This, const char* name, long* &value, int& size, int& allocated)
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
int atlas__Config__get_array_float (Config* This, const char* name, float* &value, int& size, int& allocated)
{
  ATLAS_ERROR_HANDLING(
    std::vector<double> v;
    if( ! This->get(std::string(name),v) )
      return false;
    size = v.size();
    value = new float[size];
    for ( size_t j = 0; j < v.size(); ++j ) {
      value[j] = v[j];
    }
    allocated = true;
  );
    return true;
}
int atlas__Config__get_array_double (Config* This, const char* name, double* &value, int& size, int& allocated)
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

int atlas__Config__has (Config *This, const char *name) {
    ATLAS_ERROR_HANDLING( return This->has( std::string(name) ));
    return 0;
}

void atlas__Config__json(Config* This, char* &json, int &size, int &allocated)
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

} // namespace util
} // namespace atlas
