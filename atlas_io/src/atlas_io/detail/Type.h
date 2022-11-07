/*
 * (C) Copyright 2020 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#pragma once

#include <string>

namespace atlas {
namespace io {

struct Type {
    const std::string name_;
    operator const std::string&() { return name_; }
    operator bool() const { return name_.size(); }
    Type(const char* name): name_(name) {}
    Type(const std::string& name): name_(name) {}
    bool operator==(const Type& other) const { return name_ == other.name_; }
};


}  // namespace io
}  // namespace atlas
