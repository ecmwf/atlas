#ifndef atlas_fortran_atlas_logging_h
#define atlas_fortran_atlas_logging_h
#include <stdio.h>
extern "C"
{
  int atlas__resource_int (const char* resource, int default_value);
  long atlas__resource_long (const char* resource, long default_value);
  float atlas__resource_float (const char* resource, float default_value);
  double atlas__resource_double (const char* resource, double default_value);
  const char* atlas__resource_string (const char* resource, const char* default_value);
}

#endif
