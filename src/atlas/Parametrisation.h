/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef atlas_Parametrisation_h
#define atlas_Parametrisation_h

#include <string>
#include "eckit/config/Parametrisation.h"
#include "atlas/Metadata.h"
#include "atlas/util/Debug.h"

namespace atlas {
class Grid;
class Mesh;
class FunctionSpace;

/// @brief Parametrisation class used to construct various
///        atlas components
class Parametrisation: public eckit::Parametrisation {

public:
  
// -- Constructors

  Parametrisation();

  /// @brief Constructor starting from json stream
  Parametrisation( std::istream&, const std::string &format = "json" );

  /// @brief Constructor starting from a Properties
  Parametrisation( const eckit::Properties& );

  /// @brief Constructor immediately setting a value.
  template<typename ValueT>
  Parametrisation( const std::string& name, const ValueT& value );

// -- Mutators

  /// @brief Operator that sets a key-value pair.
  /// This is useful for chaining. Together with above constructor:
  /// Parametrisation
  ///   ("key1",value1)
  ///   ("key2",value2)
  ///   ("key3",value3);
  template<typename ValueT>
  Parametrisation& operator()(const std::string& name, const ValueT& value);

  /// @brief Set a key-value parameter
  template<typename ValueT>
  Parametrisation& set(const std::string& name, const ValueT& value);
  Parametrisation& set(const std::string& name, const char* value);
  Parametrisation& set(const eckit::Properties&);
  Parametrisation& set(const Parametrisation&);

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

  bool get(const std::string& name, Parametrisation& value) const;
  bool get(const std::string& name, std::vector<Parametrisation>& value) const;

private:

  template<class T>
  bool _get(const std::string &name, T &value) const;
  Metadata delegate_;

};

// ------------------------------------------------------------------

template<typename ValueT>
Parametrisation::Parametrisation(const std::string& name, const ValueT& value)
{
  delegate_.set(name,value);
}

template<typename ValueT>
Parametrisation& Parametrisation::operator()(const std::string& name, const ValueT& value)
{
  set(name,value);
  return *this;
}

template<typename ValueT>
Parametrisation& Parametrisation::set(const std::string& name, const ValueT& value)
{
  delegate_.set(name,value);
  return *this;
}

// ------------------------------------------------------------------
// C wrapper interfaces to C++ routines

#define Char char
extern "C"
{
  Parametrisation* atlas__Parametrisation__new ();

  void atlas__Parametrisation__delete (Parametrisation* This);
  int  atlas__Parametrisation__has (Parametrisation* This, const char* name);
  void atlas__Parametrisation__set_int    (Parametrisation* This, const char* name, int value);
  void atlas__Parametrisation__set_long   (Parametrisation* This, const char* name, long value);
  void atlas__Parametrisation__set_float  (Parametrisation* This, const char* name, float value);
  void atlas__Parametrisation__set_double (Parametrisation* This, const char* name, double value);
  void atlas__Parametrisation__set_string (Parametrisation* This, const char* name, const char* value);
  void atlas__Parametrisation__set_array_int    (Parametrisation* This, const char* name, int value[], int size);
  void atlas__Parametrisation__set_array_long   (Parametrisation* This, const char* name, long value[], int size);
  void atlas__Parametrisation__set_array_float  (Parametrisation* This, const char* name, float value[], int size);
  void atlas__Parametrisation__set_array_double (Parametrisation* This, const char* name, double value[], int size);

  int atlas__Parametrisation__get_int    (Parametrisation* This, const char* name, int &value);
  int atlas__Parametrisation__get_long   (Parametrisation* This, const char* name, long &value);
  int atlas__Parametrisation__get_float  (Parametrisation* This, const char* name, float &value);
  int atlas__Parametrisation__get_double (Parametrisation* This, const char* name, double &value);
  int atlas__Parametrisation__get_string (Parametrisation* This, const char* name, Char* &value, int &size, int &allocated);
  int atlas__Parametrisation__get_array_int    (Parametrisation* This, const char* name, int* &value, int &size, int &allocated);
  int atlas__Parametrisation__get_array_long   (Parametrisation* This, const char* name, long* &value, int &size, int &allocated);
  int atlas__Parametrisation__get_array_float  (Parametrisation* This, const char* name, float* &value, int &size, int &allocated);
  int atlas__Parametrisation__get_array_double (Parametrisation* This, const char* name, double* &value, int &size, int &allocated);
}
#undef Char

// ------------------------------------------------------------------

} // namespace atlas

#endif // Parametrisation_h
