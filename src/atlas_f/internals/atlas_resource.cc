#include "eckit/config/Resource.h"
#include "eckit/config/ResourceMgr.h"
#include "eckit/utils/Translator.h"
#include "atlas_f/internals/atlas_resource.h"
#include "atlas/runtime/ErrorHandling.h"

namespace atlas {

int atlas__resource_int (const char* resource, int default_value)
{
  ATLAS_ERROR_HANDLING(
    return eckit::Resource<int>( std::string(resource), default_value )
  );
  return int(0);
}

long atlas__resource_long (const char* resource, long default_value)
{
  ATLAS_ERROR_HANDLING(
    return eckit::Resource<long>( std::string(resource), default_value )
  );
  return long(0);
}

float atlas__resource_float (const char* resource, float default_value)
{
  ATLAS_ERROR_HANDLING(
    return static_cast<float>( eckit::Resource<double>( std::string(resource), static_cast<double>(default_value) ) );
  );
  return float(0);
}

double atlas__resource_double (const char* resource, double default_value)
{
  ATLAS_ERROR_HANDLING(
    return eckit::Resource<double>( std::string(resource), default_value )
  );
  return double(0);
}

const char* atlas__resource_string (const char* resource, const char* default_value)
{
  static std::string val;
  ATLAS_ERROR_HANDLING(
    val = eckit::Resource<std::string>( std::string(resource), std::string(default_value) );
    return val.c_str();
  );
  return NULL;
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
