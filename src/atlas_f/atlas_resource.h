#ifndef atlas_fortran_atlas_logging_h
#define atlas_fortran_atlas_logging_h
#include <stdio.h>

namespace atlas{

extern "C"
{
  int atlas__resource_int (const char* resource, int default_value);
  long atlas__resource_long (const char* resource, long default_value);
  float atlas__resource_float (const char* resource, float default_value);
  double atlas__resource_double (const char* resource, double default_value);
  const char* atlas__resource_string (const char* resource, const char* default_value);
  void atlas__resource_set_int (const char* resource, int value);
  void atlas__resource_set_long (const char* resource, long value);
  void atlas__resource_set_float (const char* resource, float value);
  void atlas__resource_set_double (const char* resource, double value);
  void atlas__resource_set_string (const char* resource, const char* value);
}

}

#endif
