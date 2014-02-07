#ifndef Metadata_hpp
#define Metadata_hpp

#include <vector>
#include <map>
#include <string>

namespace ecmwf {
class Field;

/// @brief Contains a list of field-pointers, no ownership
class Metadata
{
public:
  template<typename ValueT>
  Metadata& set(const std::string& name, const ValueT& value);

  template<typename ValueT>
  ValueT& get(const std::string& name);

private:
  std::map< std::string, bool > map_bool_;
  std::map< std::string, int > map_int_;
  std::map< std::string, float > map_float_;
  std::map< std::string, double > map_double_;
  std::map< std::string, std::string > map_string_;
};

// ------------------------------------------------------------------
// C wrapper interfaces to C++ routines
extern "C" 
{
  Metadata* ecmwf__Metadata__new ();
  void ecmwf__Metadata__delete (Metadata* This);
  void ecmwf__Metadata__add_int    (Metadata* This, const char* name, int value);
  void ecmwf__Metadata__add_float  (Metadata* This, const char* name, float value);
  void ecmwf__Metadata__add_double (Metadata* This, const char* name, double value);
  void ecmwf__Metadata__add_string (Metadata* This, const char* name, const char* value);

  int    ecmwf__Metadata__get_int    (Metadata* This, const char* name);
  float  ecmwf__Metadata__get_float  (Metadata* This, const char* name);
  double ecmwf__Metadata__get_double (Metadata* This, const char* name);
  const char* ecmwf__Metadata__get_string (Metadata* This, const char* name);

}
// ------------------------------------------------------------------

} // namespace ecmwf

#endif // Metadata_hpp
