#include "eckit/config/Resource.h"
#include "atlas/fortran/atlas_resource.h"

extern "C"
{

int atlas__resource_int (const char* resource, int default_value)
{
  int val = eckit::Resource<int>( std::string(resource), default_value );
  return val;
}

long atlas__resource_long (const char* resource, long default_value)
{
  long val = eckit::Resource<long>( std::string(resource), default_value );
  return val;
}

float atlas__resource_float (const char* resource, float default_value)
{
  double val = eckit::Resource<double>( std::string(resource), static_cast<double>(default_value) );
  return static_cast<float>(val);
}

double atlas__resource_double (const char* resource, double default_value)
{
  double val = eckit::Resource<double>( std::string(resource), default_value );
  return val;
}

const char* atlas__resource_string (const char* resource, const char* default_value)
{
  static std::string val;
  val = eckit::Resource<std::string>( std::string(resource), std::string(default_value) );
  return val.c_str();
}

}
