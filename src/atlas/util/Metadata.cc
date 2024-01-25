/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/util/Metadata.h"

#include <iostream>
#include <sstream>
#include <stdexcept>

#include "eckit/log/JSON.h"
#include "eckit/parser/JSONParser.h"
#include "eckit/utils/Hash.h"


#include "atlas/parallel/mpi/mpi.h"
#include "atlas/runtime/Exception.h"

using std::string;

namespace atlas {
namespace util {

void Metadata::throw_not_found(const std::string& name) const {
    std::stringstream msg;
    msg << "Could not find metadata \"" << name << "\"";
    throw_Exception(msg.str(), Here());
}

size_t Metadata::footprint() const {
    // TODO
    size_t size = sizeof(*this);
    return size;
}

void Metadata::broadcast() {
    idx_t root = 0;
    get("owner", root);
    broadcast(*this, root);
}

void Metadata::broadcast(idx_t root) {
    broadcast(*this, root);
}

void Metadata::broadcast(Metadata& dest) {
    idx_t root = 0;
    get("owner", root);
    broadcast(dest, root);
}

void Metadata::broadcast(Metadata& dest, idx_t root) {
    std::string buffer;
    int buffer_size{0};
    if (mpi::rank() == root) {
        std::stringstream s;
        eckit::JSON json(s);
        json.precision(17);
        json << *this;
        buffer      = s.str();
        buffer_size = static_cast<int>(buffer.size());
    }

    ATLAS_TRACE_MPI(BROADCAST) { atlas::mpi::comm().broadcast(buffer_size, root); }

    if (mpi::rank() != root) {
        buffer.resize(buffer_size);
    }

    ATLAS_TRACE_MPI(BROADCAST) { mpi::comm().broadcast(buffer.begin(), buffer.end(), root); }

    if (not(&dest == this && mpi::rank() == root)) {
        std::stringstream s;
        s << buffer;
        eckit::JSONParser parser(s);
        dest = Metadata(parser.parse());
    }
}

void Metadata::broadcast(Metadata& dest) const {
    idx_t root = 0;
    get("owner", root);
    broadcast(dest, root);
}

void Metadata::broadcast(Metadata& dest, idx_t root) const {
    std::string buffer;
    int buffer_size{0};
    if (mpi::rank() == root) {
        std::stringstream s;
        eckit::JSON json(s);
        json.precision(17);
        json << *this;
        buffer      = s.str();
        buffer_size = static_cast<int>(buffer.size());
    }

    ATLAS_TRACE_MPI(BROADCAST) { mpi::comm().broadcast(buffer_size, root); }

    if (mpi::rank() != root) {
        buffer.resize(buffer_size);
    }

    ATLAS_TRACE_MPI(BROADCAST) { mpi::comm().broadcast(buffer.begin(), buffer.end(), root); }

    // Fill in dest
    {
        std::stringstream s;
        s << buffer;
        eckit::JSONParser parser(s);
        dest = Metadata(parser.parse());
    }
}


Metadata& Metadata::set(const eckit::Configuration& other) {
#if ATLAS_ECKIT_VERSION_AT_LEAST(1, 26, 0) || ATLAS_ECKIT_DEVELOP
    LocalConfiguration::set(other);
#else
    eckit::Value& root = const_cast<eckit::Value&>(get());
    auto& other_root   = other.get();
    std::vector<string> other_keys;
    eckit::fromValue(other_keys, other_root.keys());
    for (auto& key : other_keys) {
        root[key] = other_root[key];
    }
#endif
    return *this;
}

Metadata::Metadata(const eckit::Value& value): eckit::LocalConfiguration(value) {}

// ------------------------------------------------------------------
// C wrapper interfaces to C++ routines

Metadata* atlas__Metadata__new() {
    return new Metadata();
}

void atlas__Metadata__delete(Metadata* This) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_Metadata");
    delete This;
}

void atlas__Metadata__set_int(Metadata* This, const char* name, int value) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_Metadata");
    This->set(std::string(name), long(value));
}
void atlas__Metadata__set_long(Metadata* This, const char* name, long value) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_Metadata");
    This->set(std::string(name), value);
}
void atlas__Metadata__set_float(Metadata* This, const char* name, float value) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_Metadata");
    This->set(std::string(name), double(value));
}
void atlas__Metadata__set_double(Metadata* This, const char* name, double value) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_Metadata");
    This->set(std::string(name), value);
}
void atlas__Metadata__set_string(Metadata* This, const char* name, const char* value) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_Metadata");
    This->set(std::string(name), std::string(value));
}
void atlas__Metadata__set_array_int(Metadata* This, const char* name, int value[], int size) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_Metadata");
    std::vector<int> v;
    v.assign(value, value + size);
    This->set(std::string(name), v);
}
void atlas__Metadata__set_array_long(Metadata* This, const char* name, long value[], int size) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_Metadata");
    std::vector<long> v;
    v.assign(value, value + size);
    This->set(std::string(name), v);
}
void atlas__Metadata__set_array_float(Metadata* This, const char* name, float value[], int size) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_Metadata");
    std::vector<float> v;
    v.assign(value, value + size);
    This->set(std::string(name), v);
}
void atlas__Metadata__set_array_double(Metadata* This, const char* name, double value[], int size) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_Metadata");
    std::vector<double> v;
    v.assign(value, value + size);
    This->set(std::string(name), v);
}
int atlas__Metadata__get_int(Metadata* This, const char* name) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_Metadata");
    return This->get<int>(std::string(name));
}
long atlas__Metadata__get_long(Metadata* This, const char* name) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_Metadata");
    return This->get<long>(std::string(name));
}
float atlas__Metadata__get_float(Metadata* This, const char* name) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_Metadata");
    return This->get<float>(std::string(name));
}
double atlas__Metadata__get_double(Metadata* This, const char* name) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_Metadata");
    return This->get<double>(std::string(name));
}
void atlas__Metadata__get_string(Metadata* This, const char* name, char* output_str, int max_len) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_Metadata");
    std::string s = This->get<std::string>(std::string(name));
    if (s.size() + 1 > size_t(max_len)) {
        std::stringstream msg;
        msg << "Cannot copy string `" << s << "` of metadata `" << name
            << "`"
               "in buffer of length "
            << max_len;
        throw_Exception(msg.str(), Here());
    }
    std::strncpy(output_str, s.c_str(), max_len);
}
void atlas__Metadata__get_array_int(Metadata* This, const char* name, int*& value, int& size, int& allocated) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_Metadata");
    std::vector<int> v = This->get<std::vector<int>>(std::string(name));
    size               = static_cast<int>(v.size());
    value              = new int[size];
    for (size_t j = 0; j < v.size(); ++j) {
        value[j] = v[j];
    }
    allocated = true;
}
void atlas__Metadata__get_array_long(Metadata* This, const char* name, long*& value, int& size, int& allocated) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_Metadata");
    std::vector<long> v = This->get<std::vector<long>>(std::string(name));
    size                = static_cast<int>(v.size());
    value               = new long[size];
    for (size_t j = 0; j < v.size(); ++j) {
        value[j] = v[j];
    }
    allocated = true;
}
void atlas__Metadata__get_array_float(Metadata* This, const char* name, float*& value, int& size, int& allocated) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_Metadata");
    std::vector<float> v = This->get<std::vector<float>>(std::string(name));
    size                 = static_cast<int>(v.size());
    value                = new float[size];
    for (size_t j = 0; j < v.size(); ++j) {
        value[j] = v[j];
    }
    allocated = true;
}
void atlas__Metadata__get_array_double(Metadata* This, const char* name, double*& value, int& size, int& allocated) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_Metadata");
    std::vector<double> v = This->get<std::vector<double>>(std::string(name));
    size                  = static_cast<int>(v.size());
    value                 = new double[size];
    for (size_t j = 0; j < v.size(); ++j) {
        value[j] = v[j];
    }
    allocated = true;
}

int atlas__Metadata__has(Metadata* This, const char* name) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_Metadata");
    return This->has(std::string(name));
}

void atlas__Metadata__print(Metadata* This, std::ostream* channel) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialised atlas_Metadata");
    ATLAS_ASSERT(channel != nullptr);
    *channel << *This;
}

void atlas__Metadata__json(Metadata* This, char*& json, int& size, int& allocated) {
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
