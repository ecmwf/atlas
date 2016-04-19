/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef atlas_Config_h
#define atlas_Config_h

#include <string>
#include "eckit/config/Parametrisation.h"
#include "atlas/util/Metadata.h"

namespace eckit {
  class JSON;
  class PathName;
}

namespace atlas {
namespace util {
class Grid;
class Mesh;

/// @brief Parametrisation class used to construct various
///        atlas components
class Config: public eckit::Parametrisation {

public:

// -- Constructors

  Config();

  /// @brief Constructor starting from a path
  Config( const eckit::PathName& );

  /// @brief Constructor starting from a stream
  Config( std::istream&, const std::string &format = "json" );

  /// @brief Constructor starting from a Properties
  Config( const eckit::Properties& );

  /// @brief Constructor immediately setting a value.
  template<typename ValueT>
  Config( const std::string& name, const ValueT& value );

// -- Mutators

  /// @brief Operator that sets a key-value pair.
  /// This is useful for chaining. Together with above constructor:
  /// Config
  ///   ("key1",value1)
  ///   ("key2",value2)
  ///   ("key3",value3);
  template<typename ValueT>
  Config operator()(const std::string& name, const ValueT& value);
  Config operator()(const Config &c) { return set(c); }


  // Overload operators to merge two Config objects.
  Config operator&&(const Config& other) const;
  Config operator|(const Config& other) const;

  /// @brief Set a key-value parameter
  template<typename ValueT>
  Config& set(const std::string& name, const ValueT& value);
  Config& set(const std::string& name, const char* value);
  Config& set(const std::string& name, const Config& );
  Config& set(const std::string& name, const std::vector<Config>& );
  Config& set(const eckit::Properties&);
  Config& set(const Config&);

// -- Accessors, overloaded from eckit::Parametrisation

  /// @returns true if a parameter exists
  virtual bool has( const std::string& name ) const;

  virtual bool get(const std::string& name, std::string& value) const;
  virtual bool get(const std::string& name, bool& value) const;
  virtual bool get(const std::string& name, long& value) const;
  virtual bool get(const std::string& name, size_t& value) const;
  virtual bool get(const std::string& name, double& value) const;

  virtual bool get(const std::string& name, std::vector<long>& value) const;
  virtual bool get(const std::string& name, std::vector<double>& value) const;

  bool get(const std::string& name, Config& value) const;
  bool get(const std::string& name, std::vector<Config>& value) const;

  friend eckit::JSON& operator<<(eckit::JSON&, const Config&);
  friend std::ostream& operator<<(std::ostream&, const Config&);
private:

  template<class T>
  bool _get(const std::string &name, T &value) const;
  Metadata delegate_;

};

// ------------------------------------------------------------------

template<typename ValueT>
inline Config::Config(const std::string& name, const ValueT& value)
{
  delegate_.set(name,value);
}

template<typename ValueT>
inline Config Config::operator()(const std::string& name, const ValueT& value)
{
  set(name,value);
  return *this;
}

template<typename ValueT>
inline Config& Config::set(const std::string& name, const ValueT& value)
{
  delegate_.set(name,value);
  return *this;
}


// ------------------------------------------------------------------
// C wrapper interfaces to C++ routines

#define Char char
extern "C"
{
  Config* atlas__Config__new ();
  Config* atlas__Config__new_from_json (const char* json);
  Config* atlas__Config__new_from_file (const char* path);
  void atlas__Config__delete (Config* This);
  int  atlas__Config__has (Config* This, const char* name);
  void atlas__Config__set_config (Config* This, const char* name, const Config* value);
  void atlas__Config__set_config_list (Config* This, const char* name, const Config* value[], int size);
  void atlas__Config__set_int    (Config* This, const char* name, int value);
  void atlas__Config__set_long   (Config* This, const char* name, long value);
  void atlas__Config__set_float  (Config* This, const char* name, float value);
  void atlas__Config__set_double (Config* This, const char* name, double value);
  void atlas__Config__set_string (Config* This, const char* name, const char* value);
  void atlas__Config__set_array_int    (Config* This, const char* name, int value[], int size);
  void atlas__Config__set_array_long   (Config* This, const char* name, long value[], int size);
  void atlas__Config__set_array_float  (Config* This, const char* name, float value[], int size);
  void atlas__Config__set_array_double (Config* This, const char* name, double value[], int size);

  int atlas__Config__get_config (Config* This, const char* name, Config* value);
  int atlas__Config__get_config_list (Config* This, const char* name, Config** &value, int &size, int &allocated);
  int atlas__Config__get_int    (Config* This, const char* name, int &value);
  int atlas__Config__get_long   (Config* This, const char* name, long &value);
  int atlas__Config__get_float  (Config* This, const char* name, float &value);
  int atlas__Config__get_double (Config* This, const char* name, double &value);
  int atlas__Config__get_string (Config* This, const char* name, Char* &value, int &size, int &allocated);
  int atlas__Config__get_array_int    (Config* This, const char* name, int* &value, int &size, int &allocated);
  int atlas__Config__get_array_long   (Config* This, const char* name, long* &value, int &size, int &allocated);
  int atlas__Config__get_array_float  (Config* This, const char* name, float* &value, int &size, int &allocated);
  int atlas__Config__get_array_double (Config* This, const char* name, double* &value, int &size, int &allocated);
  void atlas__Config__json (Config* This, Char* &json, int &size, int &allocated);
}
#undef Char

// ------------------------------------------------------------------

class NoConfig: public eckit::Parametrisation {

public: // methods

    /// Destructor redundant but fixes sanity compiler warnings
    virtual ~NoConfig() {}

    virtual bool has(const std::string& name) const { return false; }

    virtual bool get(const std::string& name, std::string& value) const { return false; }
    virtual bool get(const std::string& name, bool& value) const { return false; }
    virtual bool get(const std::string& name, long& value) const { return false; }
    virtual bool get(const std::string& name, size_t& value) const { return false; }
    virtual bool get(const std::string& name, double& value) const { return false; }

    virtual bool get(const std::string& name, std::vector<long>& value) const { return false; }
    virtual bool get(const std::string& name, std::vector<double>& value) const { return false; }

};

} // namespace util
} // namespace atlas

#endif // Parametrisation_h
