/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "atlas/util/Metadata.h"

#include <iostream>
#include <stdexcept>
#include <sstream>

#include "eckit/exception/Exceptions.h"
#include "atlas/parallel/mpi/mpi.h"
#include "eckit/parser/JSON.h"
#include "eckit/parser/JSONParser.h"
#include "eckit/utils/Hash.h"

#include "atlas/runtime/ErrorHandling.h"

using std::string;

namespace atlas {
namespace util {

void Metadata::throw_exception(const std::string& name) const {
      std::stringstream msg;
      msg << "Could not find metadata \"" << name << "\"";
      throw eckit::OutOfRange(msg.str(),Here());
}

size_t Metadata::footprint() const
{
  // TODO
  size_t size = sizeof(*this);
  return size;
}

void Metadata::broadcast()
{
  size_t root = 0;
  get( "owner", root );
  broadcast(*this,root);
}

void Metadata::broadcast(const size_t root)
{
  broadcast(*this,root);
}

void Metadata::broadcast(Metadata& dest)
{
  size_t root = 0;
  get( "owner", root );
  broadcast(dest,root);
}

void Metadata::broadcast(Metadata& dest, const size_t root)
{
  std::string buffer;
  int buffer_size;
  if( atlas::parallel::mpi::comm().rank() == root )
  {
    std::stringstream s;
    eckit::JSON json(s);
    json.precision(17);
    json << *this;
    buffer = s.str();
    buffer_size = buffer.size();
  }

  ATLAS_TRACE_MPI( BROADCAST ) {
    atlas::parallel::mpi::comm().broadcast(buffer_size,root);
  }

  if( atlas::parallel::mpi::comm().rank() != root ) {
    buffer.resize(buffer_size);
  }

  ATLAS_TRACE_MPI( BROADCAST ) {
    atlas::parallel::mpi::comm().broadcast(buffer.begin(), buffer.end(), root);
  }

  if( not (&dest==this && atlas::parallel::mpi::comm().rank() == root ) )
  {
    std::stringstream s;
    s << buffer;
    eckit::JSONParser parser( s );
    dest = Metadata( parser.parse() );
  }
}

void Metadata::broadcast(Metadata& dest) const
{
  size_t root = 0;
  get( "owner", root );
  broadcast(dest,root);
}

void Metadata::broadcast(Metadata& dest, const size_t root) const
{
  std::string buffer;
  int buffer_size;
  if( atlas::parallel::mpi::comm().rank() == root )
  {
    std::stringstream s;
    eckit::JSON json(s);
    json.precision(17);
    json << *this;
    buffer = s.str();
    buffer_size = buffer.size();
  }

  ATLAS_TRACE_MPI( BROADCAST ) {
    atlas::parallel::mpi::comm().broadcast(buffer_size,root);
  }

  if( atlas::parallel::mpi::comm().rank() != root ) {
    buffer.resize(buffer_size);
  }

  ATLAS_TRACE_MPI( BROADCAST ) {
    atlas::parallel::mpi::comm().broadcast(buffer.begin(), buffer.end(), root);
  }

  // Fill in dest
  {
    std::stringstream s;
    s << buffer;
    eckit::JSONParser parser( s );
    dest = Metadata( parser.parse() );
  }
}


void Metadata::hash(eckit::Hash& hsh) const
{
  eckit::ValueMap map = get();
  for( eckit::ValueMap::const_iterator vit = map.begin(); vit != map.end(); ++vit )
  {
    hsh.add(vit->first.as<std::string>());
    /// @note below, we assume all Values translate to std::string, this needs more verification
    hsh.add(vit->second.as<std::string>());
  }
}

Metadata::Metadata( const eckit::Value& value ) :
    eckit::LocalConfiguration(value) {
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

} // namespace util
} // namespace atlas
