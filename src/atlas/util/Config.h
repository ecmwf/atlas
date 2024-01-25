/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#pragma once
#include <initializer_list>
#include <string>
#include <string_view>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/log/JSON.h"

#include "atlas/library/config.h"

namespace eckit {
class PathName;
class Hash;
}  // namespace eckit

namespace atlas {
namespace util {

/// @brief Configuration class used to construct various
///        atlas components
class Config : public eckit::LocalConfiguration {
public:
    // -- Constructors

    Config();

    /// @brief Constructor starting from a path
    Config(const eckit::PathName&);

    /// @brief Constructor starting from a stream
    Config(std::istream&, const std::string& format = "json");

    /// @brief Constructor starting from a Configuration
    Config(const eckit::Configuration&);

    /// @brief Constructor immediately setting a value.
    template <typename ValueT>
    Config(const std::string& name, const ValueT& value);

    Config(const std::string& name, const std::string_view value) :
        Config(name, std::string(value)) {}

    template <typename ValueT>
    Config(const std::string& name, std::initializer_list<ValueT>&& value);

    // -- Mutators

    /// @brief Operator that sets a key-value pair.
    template <typename ValueT>
    Config operator()(const std::string& name, const ValueT& value);

    template <typename ValueT>
    Config operator()(const std::string& name, std::initializer_list<ValueT>&& value);

    // Overload operators to merge two Config objects.
    Config operator|(const eckit::Configuration& other) const;

    /// @brief Set a key-value parameter
    using eckit::LocalConfiguration::set;

    Config& set(const std::string& name, const std::vector<Config>&);

    Config& set(const eckit::Configuration&);

    template <typename T>
    Config& set(const std::string& name, std::initializer_list<T>&& value);

    Config& remove(const std::string& name);

    // -- Accessors, overloaded from eckit::Parametrisation

    using eckit::LocalConfiguration::get;
    bool get(const std::string& name, std::vector<Config>& value) const;

#if ! (ATLAS_ECKIT_VERSION_AT_LEAST(1, 26, 0) || ATLAS_ECKIT_DEVELOP)
    std::vector<std::string> keys() const;
#endif

    std::string json(eckit::JSON::Formatting = eckit::JSON::Formatting::indent()) const;
};

// ------------------------------------------------------------------

template <typename ValueT>
inline Config::Config(const std::string& name, const ValueT& value): eckit::LocalConfiguration() {
    set(name, value);
}

template <typename ValueT>
inline Config::Config(const std::string& name, std::initializer_list<ValueT>&& value): eckit::LocalConfiguration() {
    set(name, value);
}

template <typename ValueT>
inline Config Config::operator()(const std::string& name, const ValueT& value) {
    set(name, value);
    return *this;
}

template <typename ValueT>
inline Config& Config::set(const std::string& name, std::initializer_list<ValueT>&& value) {
    set(name, std::vector<ValueT>(value.begin(), value.end()));
    return *this;
}

template <typename ValueT>
inline Config Config::operator()(const std::string& name, std::initializer_list<ValueT>&& value) {
    set(name, std::vector<ValueT>(value.begin(), value.end()));
    return *this;
}

// ------------------------------------------------------------------
// C wrapper interfaces to C++ routines
// clang-format off
extern "C" {
Config* atlas__Config__new();
Config* atlas__Config__new_from_json( const char* json );
Config* atlas__Config__new_from_file( const char* path );
void atlas__Config__delete( Config* This );
int atlas__Config__has( Config* This, const char* name );
void atlas__Config__set_config( Config* This, const char* name, const Config* value );
void atlas__Config__set_config_list( Config* This, const char* name, const Config* value[], int size );
void atlas__Config__set_int( Config* This, const char* name, int value );
void atlas__Config__set_long( Config* This, const char* name, long value );
void atlas__Config__set_float( Config* This, const char* name, float value );
void atlas__Config__set_double( Config* This, const char* name, double value );
void atlas__Config__set_string( Config* This, const char* name, const char* value );
void atlas__Config__set_array_int( Config* This, const char* name, int value[], int size );
void atlas__Config__set_array_long( Config* This, const char* name, long value[], int size );
void atlas__Config__set_array_float( Config* This, const char* name, float value[], int size );
void atlas__Config__set_array_double( Config* This, const char* name, double value[], int size );

int atlas__Config__get_config( Config* This, const char* name, Config* value );
int atlas__Config__get_config_list( Config* This, const char* name, Config** &value, int& size, int& allocated );
int atlas__Config__get_int( Config* This, const char* name, int& value );
int atlas__Config__get_long( Config* This, const char* name, long& value );
int atlas__Config__get_float( Config* This, const char* name, float& value );
int atlas__Config__get_double( Config* This, const char* name, double& value );
int atlas__Config__get_string( Config* This, const char* name, char*& value, int& size, int& allocated );
int atlas__Config__get_array_int( Config* This, const char* name, int*& value, int& size, int& allocated );
int atlas__Config__get_array_long( Config* This, const char* name, long*& value, int& size, int& allocated );
int atlas__Config__get_array_float( Config* This, const char* name, float*& value, int& size, int& allocated );
int atlas__Config__get_array_double( Config* This, const char* name, double*& value, int& size, int& allocated );
void atlas__Config__json( Config* This, char*& json, int& size, int& allocated );
}
// clang-format on
// ------------------------------------------------------------------

class NoConfig : public Config {
public:  // methods
    /// Destructor redundant but fixes sanity compiler warnings
    virtual ~NoConfig() {}

    virtual bool has(const std::string& name) const { return false; }

    virtual bool get(const std::string& name, std::string& value) const { return false; }
    virtual bool get(const std::string& name, bool& value) const { return false; }
    virtual bool get(const std::string& name, int& value) const { return false; }
    virtual bool get(const std::string& name, long& value) const { return false; }
    virtual bool get(const std::string& name, long long& value) const { return false; }
    virtual bool get(const std::string& name, std::size_t& value) const { return false; }
    virtual bool get(const std::string& name, float& value) const { return false; }
    virtual bool get(const std::string& name, double& value) const { return false; }

    virtual bool get(const std::string& name, std::vector<std::string>& value) const { return false; }
    virtual bool get(const std::string& name, std::vector<int>& value) const { return false; }
    virtual bool get(const std::string& name, std::vector<long>& value) const { return false; }
    virtual bool get(const std::string& name, std::vector<long long>& value) const { return false; }
    virtual bool get(const std::string& name, std::vector<std::size_t>& value) const { return false; }
    virtual bool get(const std::string& name, std::vector<float>& value) const { return false; }
    virtual bool get(const std::string& name, std::vector<double>& value) const { return false; }
};

}  // namespace util
}  // namespace atlas
