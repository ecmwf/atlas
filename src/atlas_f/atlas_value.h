#ifndef atlas_f_atlas_value_h
#define atlas_f_atlas_value_h

namespace eckit{ class Value; }
#define AtlasValue eckit::Value
#define CHAR char
namespace atlas{

extern "C"
{
  AtlasValue* atlas__Value__new_int (int value);
  AtlasValue* atlas__Value__new_long (long value);
  AtlasValue* atlas__Value__new_float (float value);
  AtlasValue* atlas__Value__new_double (double value);
  AtlasValue* atlas__Value__new_string (const char* value);
  AtlasValue* atlas__Value__new_array_int (int value[], int size);
  AtlasValue* atlas__Value__new_array_long (long value[], int size);
  AtlasValue* atlas__Value__new_array_float (float value[], int size);
  AtlasValue* atlas__Value__new_array_double (double value[], int size);
  void atlas__Value__delete (AtlasValue* This);

  void atlas__Value__int (AtlasValue* This, int &value);
  void atlas__Value__long (AtlasValue* This, long &value);
  void atlas__Value__float (AtlasValue* This, float &value);
  void atlas__Value__double (AtlasValue* This, double &value);
  void atlas__Value__string (AtlasValue* This, CHAR* &value, int &size, int &allocated);
  void atlas__Value__array_int (AtlasValue* This, int* &value, int &size, int &allocated);
  void atlas__Value__array_long (AtlasValue* This, long* &value, int &size, int &allocated);
  void atlas__Value__array_float (AtlasValue* This, float* &value, int &size, int &allocated);
  void atlas__Value__array_double (AtlasValue* This, double* &value, int &size, int &allocated);
}

}
#undef AtlasValue
#undef CHAR
#endif
