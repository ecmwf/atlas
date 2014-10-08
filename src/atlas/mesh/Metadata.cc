/*
 * (C) Copyright 1996-2014 ECMWF.
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
#include "eckit/exception/Exceptions.h"
#include "atlas/mesh/Metadata.h"
#include "atlas/mesh/Field.h"

using namespace std;

#define METADATA( VALUE_TYPE ) \
template<>\
bool Metadata::has<VALUE_TYPE>(const std::string& name) const\
{\
	return map_##VALUE_TYPE##_.count(name);\
}\
template<>\
Metadata& Metadata::set(const std::string& name, const VALUE_TYPE& value)\
{\
	map_##VALUE_TYPE##_[name] = value;\
	return *this;\
}\
template<>\
const VALUE_TYPE& Metadata::get(const std::string& name) const\
{\
	if( has<VALUE_TYPE>(name) ) {\
		return map_##VALUE_TYPE##_.at(name);\
	}\
	else {\
		std::stringstream msg;\
		msg << "Could not find metadata \"" << name << "\"";\
		throw eckit::OutOfRange(msg.str(),Here());\
	}\
}\
template<>\
const VALUE_TYPE& Metadata::get(const std::string& name, const VALUE_TYPE& deflault_value) const\
{\
	if( has<VALUE_TYPE>(name) ) {\
		return map_##VALUE_TYPE##_.at(name);\
	}\
	else {\
		return deflault_value;\
	}\
}


#define METADATA_C_BINDING( VALUE_TYPE ) \
void atlas__Metadata__add_##VALUE_TYPE (Metadata* This, const char* name, VALUE_TYPE value)\
{\
	This->set( std::string(name) ,value);\
}\
VALUE_TYPE atlas__Metadata__get_##VALUE_TYPE (Metadata* This, const char* name)\
{\
	return This->get<VALUE_TYPE>( std::string(name) );\
}


namespace atlas {

METADATA(bool)
METADATA(int)
METADATA(float)
METADATA(double)
METADATA(string)
//METADATA(void*)

// ------------------------------------------------------------------
// C wrapper interfaces to C++ routines

Metadata* atlas__Metadata__new () {
	return new Metadata();
}

void atlas__Metadata__delete (Metadata* This) {
	delete This;
}

METADATA_C_BINDING(int)
METADATA_C_BINDING(float)
METADATA_C_BINDING(double)

void atlas__Metadata__add_string (Metadata* This, const char* name, const char* value)
{
	This->set( std::string(name), std::string(value) );
}
const char* atlas__Metadata__get_string (Metadata* This, const char* name)
{
	return This->get<std::string>( std::string(name) ).c_str();
}
int	atlas__Metadata__has (Metadata* This, const char* name)
{
	return (
		This->has<int>( std::string(name) ) ||
		This->has<float>( std::string(name) ) ||
		This->has<double>( std::string(name) ) ||
		This->has<std::string>( std::string(name) ) );
}

// ------------------------------------------------------------------


} // namespace atlas

#undef METADATA
#undef METADATA_C_BINDING


