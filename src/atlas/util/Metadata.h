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

#include <string>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/config/Parametrisation.h"

#include "atlas/library/config.h"

namespace atlas {
namespace util {

class Metadata : public eckit::LocalConfiguration {
public:
    Metadata(): eckit::LocalConfiguration() {}

    Metadata(const eckit::LocalConfiguration& p): eckit::LocalConfiguration(p) {}
    //Metadata(const Metadata& p): eckit::LocalConfiguration(p) {}

    using eckit::LocalConfiguration::set;

    template <typename ValueT>
    Metadata& set(const std::string& name, const ValueT& value) {
        eckit::LocalConfiguration::set(name, value);
        return *this;
    }

    Metadata& set(const eckit::Configuration&);

    using eckit::LocalConfiguration::get;

    template <typename ValueT>
    ValueT get(const std::string& name) const {
        ValueT value;
        if (not eckit::LocalConfiguration::get(name, value))
            throw_not_found(name);
        return value;
    }

    friend std::ostream& operator<<(std::ostream& s, const Metadata& v) {
        v.print(s);
        return s;
    }

    void broadcast();
    void broadcast(idx_t root);
    void broadcast(Metadata&);
    void broadcast(Metadata&, idx_t root);
    void broadcast(Metadata&) const;
    void broadcast(Metadata&, idx_t root) const;

    size_t footprint() const;

private:
    [[noreturn]] void throw_not_found(const std::string&) const;

    Metadata(const eckit::Value&);
};

// ------------------------------------------------------------------
// C wrapper interfaces to C++ routines

extern "C" {
Metadata* atlas__Metadata__new();

void atlas__Metadata__delete(Metadata* This);
void atlas__Metadata__print(Metadata* This, std::ostream* channel);
void atlas__Metadata__json(Metadata* This, char*& json, int& size, int& allocated);
int atlas__Metadata__has(Metadata* This, const char* name);
void atlas__Metadata__set_int(Metadata* This, const char* name, int value);
void atlas__Metadata__set_long(Metadata* This, const char* name, long value);
void atlas__Metadata__set_float(Metadata* This, const char* name, float value);
void atlas__Metadata__set_double(Metadata* This, const char* name, double value);
void atlas__Metadata__set_string(Metadata* This, const char* name, const char* value);
void atlas__Metadata__set_array_int(Metadata* This, const char* name, int value[], int size);
void atlas__Metadata__set_array_long(Metadata* This, const char* name, long value[], int size);
void atlas__Metadata__set_array_float(Metadata* This, const char* name, float value[], int size);
void atlas__Metadata__set_array_double(Metadata* This, const char* name, double value[], int size);

int atlas__Metadata__get_int(Metadata* This, const char* name);
long atlas__Metadata__get_long(Metadata* This, const char* name);
float atlas__Metadata__get_float(Metadata* This, const char* name);
double atlas__Metadata__get_double(Metadata* This, const char* name);
void atlas__Metadata__get_string(Metadata* This, const char* name, char* output_str, int max_len);
void atlas__Metadata__get_array_int(Metadata* This, const char* name, int*& value, int& size, int& allocated);
void atlas__Metadata__get_array_long(Metadata* This, const char* name, long*& value, int& size, int& allocated);
void atlas__Metadata__get_array_float(Metadata* This, const char* name, float*& value, int& size, int& allocated);
void atlas__Metadata__get_array_double(Metadata* This, const char* name, double*& value, int& size, int& allocated);
}

// ------------------------------------------------------------------

}  // namespace util
}  // namespace atlas
