#include "Metadata.hpp"
#include "Field.hpp"
using namespace std;

#define METADATA( VALUE_TYPE ) \
template<>\
Metadata& Metadata::add(const std::string& name, const VALUE_TYPE& value)\
{\
  map_##VALUE_TYPE##_[name] = value;\
  return *this;\
}\
template<>\
VALUE_TYPE& Metadata::get(const std::string& name)\
{\
  return map_##VALUE_TYPE##_.at(name);\
}

#define METADATA_C_BINDING( VALUE_TYPE ) \
void ecmwf__Metadata__add_##VALUE_TYPE (Metadata* This, const char* name, VALUE_TYPE value)\
{\
  This->add( std::string(name) ,value);\
}\
VALUE_TYPE ecmwf__Metadata__get_##VALUE_TYPE (Metadata* This, const char* name)\
{\
  return This->get<VALUE_TYPE>( std::string(name) );\
}


namespace ecmwf {

METADATA(bool)
METADATA(int)
METADATA(float)
METADATA(double)
METADATA(string)

// ------------------------------------------------------------------
// C wrapper interfaces to C++ routines

Metadata* ecmwf__Metadata__new () { 
  return new Metadata(); 
}

void ecmwf__Metadata__delete (Metadata* This) {
  delete This;
}

METADATA_C_BINDING(int)
METADATA_C_BINDING(float)
METADATA_C_BINDING(double)

void ecmwf__Metadata__add_string (Metadata* This, const char* name, const char* value)
{
  This->add( std::string(name), std::string(value) );
}
const char* ecmwf__Metadata__get_string (Metadata* This, const char* name)
{
  return This->get<std::string>( std::string(name) ).c_str();
}
// ------------------------------------------------------------------


} // namespace ecmwf

#undef METADATA
#undef METADATA_C_BINDING


