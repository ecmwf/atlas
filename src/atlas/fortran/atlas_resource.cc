#include "eckit/config/Resource.h"
#include "eckit/config/ResourceMgr.h"
#include "eckit/utils/Translator.h"
#include "atlas/fortran/atlas_resource.h"
#include "atlas/ErrorHandling.h"

namespace atlas {

int atlas__resource_int (const char* resource, int default_value)
{
  int val;
  ATLAS_ERROR_HANDLING(
    val = eckit::Resource<int>( std::string(resource), default_value )
  );
  return val;
}

long atlas__resource_long (const char* resource, long default_value)
{
  long val;
  ATLAS_ERROR_HANDLING(
    val = eckit::Resource<long>( std::string(resource), default_value )
  );
  return val;
}

float atlas__resource_float (const char* resource, float default_value)
{
  double val;
  ATLAS_ERROR_HANDLING(
    val = eckit::Resource<double>( std::string(resource), static_cast<double>(default_value) )
  );
  return static_cast<float>(val);
}

double atlas__resource_double (const char* resource, double default_value)
{
  double val;
  ATLAS_ERROR_HANDLING(
    val = eckit::Resource<double>( std::string(resource), default_value )
  );
  return val;
}

const char* atlas__resource_string (const char* resource, const char* default_value)
{
  static std::string val;
  ATLAS_ERROR_HANDLING(
    val = eckit::Resource<std::string>( std::string(resource), std::string(default_value) )
  );
  return val.c_str();
}

void atlas__resource_set_int (const char* resource, int value)
{
  ATLAS_ERROR_HANDLING(
    std::stringstream val_str; val_str << value;
    eckit::ResourceMgr::instance().set(resource,val_str.str());
  );
}

void atlas__resource_set_long (const char* resource, long value)
{
  ATLAS_ERROR_HANDLING(
    std::stringstream val_str; val_str << value;
    eckit::ResourceMgr::instance().set(resource,val_str.str());
  );
}

void atlas__resource_set_float (const char* resource, float value)
{
  ATLAS_ERROR_HANDLING(
    std::stringstream val_str; val_str << value;
    eckit::ResourceMgr::instance().set(resource,val_str.str());
  );
}

void atlas__resource_set_double (const char* resource, double value)
{
  ATLAS_ERROR_HANDLING(
    std::stringstream val_str; val_str << value;
    eckit::ResourceMgr::instance().set(resource,val_str.str());
  );
}

void atlas__resource_set_string (const char* resource, const char* value)
{
  ATLAS_ERROR_HANDLING(
    std::stringstream val_str; val_str << value;
    eckit::ResourceMgr::instance().set(resource,val_str.str());
  );
}

}
