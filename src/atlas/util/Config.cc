/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <cstdarg>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include "eckit/filesystem/PathName.h"
#include "eckit/log/JSON.h"
#include "eckit/parser/YAMLParser.h"

#include "atlas/grid/Grid.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/runtime/Exception.h"
#include "atlas/util/Config.h"

using std::string;

namespace atlas {
namespace util {

namespace {

static eckit::Value yaml_from_stream(std::istream& stream) {
    eckit::YAMLParser parser(stream);
    return parser.parse();
}

static eckit::Value yaml_from_path(const eckit::PathName& path) {
    if (!path.exists()) {
        throw_Exception("File " + std::string(path) + " does not exist.");
    }
    std::ifstream file(path.localPath());
    if (!file.is_open()) {
        throw_Exception("Unable to open json file " + std::string(path), Here());
    }
    eckit::Value value = yaml_from_stream(file);
    file.close();
    return value;
}

}  // namespace

Config::Config(): eckit::LocalConfiguration() {}

Config::Config(const eckit::Configuration& p): eckit::LocalConfiguration(p) {}

Config::Config(std::istream& stream, const std::string&): eckit::LocalConfiguration(yaml_from_stream(stream)) {}

Config::Config(const eckit::PathName& path): eckit::LocalConfiguration(yaml_from_path(path)) {}

Config Config::operator|(const Config& other) const {
    Config config(*this);
    config.set(other);
    return config;
}

Config& Config::set(const eckit::LocalConfiguration& other) {
    eckit::Value& root = const_cast<eckit::Value&>(get());
    auto& other_root   = other.get();
    std::vector<string> other_keys;
    eckit::fromValue(other_keys, other_root.keys());
    for (auto& key : other_keys) {
        root[key] = other_root[key];
    }
    return *this;
}

Config& Config::remove(const std::string& name) {
    eckit::Value& root = const_cast<eckit::Value&>(get());
    root.remove(name);
    return *this;
}

Config& Config::set(const std::string& name, const std::vector<Config>& values) {
    std::vector<eckit::LocalConfiguration> metadatavalues(values.size());
    for (size_t i = 0; i < metadatavalues.size(); ++i) {
        metadatavalues[i] = values[i];
    }
    set(name, metadatavalues);
    return *this;
}

bool Config::get(const std::string& name, std::vector<Config>& value) const {
    bool found = has(name);
    if (found) {
        std::vector<eckit::LocalConfiguration> properties = getSubConfigurations(name);
        value.resize(properties.size());
        for (size_t i = 0; i < value.size(); ++i) {
            value[i] = Config(properties[i]);
        }
    }
    return found;
}

std::vector<std::string> Config::keys() const {
    std::vector<std::string> result;
    eckit::fromValue(result, get().keys());
    return result;
}

std::string Config::json(eckit::JSON::Formatting formatting) const {
    std::stringstream json;
    eckit::JSON js(json,formatting);
    js << *this;
    return json.str();
}

//==================================================================

// ------------------------------------------------------------------
// C wrapper interfaces to C++ routines

Config* atlas__Config__new() {
    return new Config();
}

Config* atlas__Config__new_from_json(const char* json) {
    std::stringstream s;
    s << json;
    return new Config(s);
}

Config* atlas__Config__new_from_file(const char* path) {
    return new Config(eckit::PathName(path));
}

void atlas__Config__delete(Config* This) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_Config");
    delete This;
}

void atlas__Config__set_config(Config* This, const char* name, const Config* value) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_Config");
    This->set(std::string(name), *value);
}

void atlas__Config__set_config_list(Config* This, const char* name, const Config* value[], int size) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_Config");
    std::vector<Config> params(size);
    for (int i = 0; i < size; ++i) {
        params[i] = Config(*value[i]);
    }
    This->set(std::string(name), params);
}

void atlas__Config__set_int(Config* This, const char* name, int value) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_Config");
    This->set(std::string(name), long(value));
}
void atlas__Config__set_long(Config* This, const char* name, long value) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_Config");
    This->set(std::string(name), value);
}
void atlas__Config__set_float(Config* This, const char* name, float value) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_Config");
    This->set(std::string(name), double(value));
}
void atlas__Config__set_double(Config* This, const char* name, double value) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_Config");
    This->set(std::string(name), value);
}
void atlas__Config__set_string(Config* This, const char* name, const char* value) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_Config");
    This->set(std::string(name), std::string(value));
}
void atlas__Config__set_array_int(Config* This, const char* name, int value[], int size) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_Config");
    std::vector<int> v;
    v.assign(value, value + size);
    This->set(std::string(name), v);
}
void atlas__Config__set_array_long(Config* This, const char* name, long value[], int size) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_Config");
    std::vector<long> v;
    v.assign(value, value + size);
    This->set(std::string(name), v);
}
void atlas__Config__set_array_float(Config* This, const char* name, float value[], int size) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_Config");
    std::vector<float> v;
    v.assign(value, value + size);
    This->set(std::string(name), v);
}
void atlas__Config__set_array_double(Config* This, const char* name, double value[], int size) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_Config");
    std::vector<double> v;
    v.assign(value, value + size);
    This->set(std::string(name), v);
}

int atlas__Config__get_config(Config* This, const char* name, Config* value) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_Config");
    if (!This->get(std::string(name), *value)) {
        return false;
    }
    return true;
}

int atlas__Config__get_config_list(Config* This, const char* name, Config**& value, int& size, int& allocated) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_Config");
    value = nullptr;
    std::vector<Config> vector;
    if (!This->get(std::string(name), vector)) {
        return false;
    }
    size      = vector.size();
    value     = new Config*[size];
    allocated = true;
    for (int i = 0; i < size; ++i) {
        value[i] = new Config(vector[i]);
    }
    return true;
}

int atlas__Config__get_int(Config* This, const char* name, int& value) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_Config");
    long long_value = value;
    if (!This->get(std::string(name), long_value)) {
        return false;
    }
    ATLAS_ASSERT(int(long_value) == long_value);
    value = long_value;
    return true;
}
int atlas__Config__get_long(Config* This, const char* name, long& value) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_Config");
    if (!This->get(std::string(name), value)) {
        return false;
    }
    return true;
}
int atlas__Config__get_float(Config* This, const char* name, float& value) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_Config");
    double double_value;
    if (!This->get(std::string(name), double_value)) {
        return false;
    }
    value = static_cast<float>(double_value);
    return true;
}
int atlas__Config__get_double(Config* This, const char* name, double& value) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_Config");
    if (!This->get(std::string(name), value)) {
        return false;
    }
    return true;
}
int atlas__Config__get_string(Config* This, const char* name, char*& value, int& size, int& allocated) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_Config");
    std::string s;
    if (!This->get(std::string(name), s)) {
        value = nullptr;
        return false;
    }
    size  = static_cast<int>(s.size());
    value = new char[size + 1];
    std::strncpy(value, s.c_str(), size + 1);
    allocated = true;
    return true;
}
int atlas__Config__get_array_int(Config* This, const char* name, int*& value, int& size, int& allocated) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_Config");
    std::vector<long> v;
    if (!This->get(std::string(name), v)) {
        return false;
    }
    size  = v.size();
    value = new int[size];
    for (size_t j = 0; j < v.size(); ++j) {
        ATLAS_ASSERT(int(v[j]) == v[j]);
        value[j] = v[j];
    }
    allocated = true;
    return true;
}
int atlas__Config__get_array_long(Config* This, const char* name, long*& value, int& size, int& allocated) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_Config");
    std::vector<long> v;
    if (!This->get(std::string(name), v)) {
        return false;
    }
    size  = v.size();
    value = new long[size];
    for (size_t j = 0; j < v.size(); ++j) {
        value[j] = v[j];
    }
    allocated = true;
    return true;
}
int atlas__Config__get_array_float(Config* This, const char* name, float*& value, int& size, int& allocated) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_Config");
    std::vector<double> v;
    if (!This->get(std::string(name), v)) {
        return false;
    }
    size  = v.size();
    value = new float[size];
    for (size_t j = 0; j < v.size(); ++j) {
        value[j] = static_cast<float>(v[j]);
    }
    allocated = true;
    return true;
}
int atlas__Config__get_array_double(Config* This, const char* name, double*& value, int& size, int& allocated) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_Config");
    std::vector<double> v;
    if (!This->get(std::string(name), v)) {
        return false;
    }
    size  = v.size();
    value = new double[size];
    for (size_t j = 0; j < v.size(); ++j) {
        value[j] = v[j];
    }
    allocated = true;
    return true;
}

int atlas__Config__has(Config* This, const char* name) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_Config");
    return This->has(std::string(name));
}

void atlas__Config__json(Config* This, char*& json, int& size, int& allocated) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_Config");
    std::stringstream s;
    eckit::JSON j(s);
    j.precision(16);
    j << *This;
    std::string json_str = s.str();
    size                 = static_cast<int>(json_str.size());
    json                 = new char[size + 1];
    allocated            = true;
    std::strncpy(json, json_str.c_str(), size + 1);
    allocated = true;
}

// ------------------------------------------------------------------

}  // namespace util
}  // namespace atlas
