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

#include "atlas/util/DataType.h"

#include "atlas/util/Config.h"

// ----------------------------------------------------------------------------

namespace atlas {
namespace option {

// ----------------------------------------------------------------------------

class type : public util::Config {
public:
    type(const std::string&);
};

// ----------------------------------------------------------------------------

class global : public util::Config {
public:
    global(size_t owner = 0);
};

// ----------------------------------------------------------------------------

class levels : public util::Config {
public:
    levels(size_t);
};

// ----------------------------------------------------------------------------

class variables : public util::Config {
public:
    variables(size_t);
};

// ----------------------------------------------------------------------------

class vector : public util::Config {
public:
    vector(size_t = 2);
};

// ----------------------------------------------------------------------------

class name : public util::Config {
public:
    name(const std::string&);
};

// ----------------------------------------------------------------------------

template <typename T>
class datatypeT : public util::Config {
public:
    datatypeT();
};

// ----------------------------------------------------------------------------

class datatype : public util::Config {
public:
    datatype(DataType::kind_t);
    datatype(const std::string&);
    datatype(DataType);
};

// ----------------------------------------------------------------------------

class shape : public util::Config {
public:
    template <typename T>
    shape(const std::initializer_list<T>& list) {
        set("shape", list);
    }
};

class alignment : public util::Config {
public:
    alignment(int);
};

// ----------------------------------------------------------------------------

class halo : public util::Config {
public:
    halo(size_t size);
};

// ----------------------------------------------------------------------------

class radius : public util::Config {
public:
    radius(double);
    radius(const std::string& = "Earth");
};

// ---------------------------------

class pole_edges : public util::Config {
public:
    pole_edges(bool = true);
};

// ----------------------------------------------------------------------------
// Definitions
// ----------------------------------------------------------------------------

template <typename T>
datatypeT<T>::datatypeT() {
    set("datatype", DataType::kind<T>());
}

}  // namespace option
}  // namespace atlas
