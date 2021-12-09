/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "eckit/exception/Exceptions.h"

#include "atlas/runtime/Exception.h"

namespace atlas {

void throw_NotImplemented(const eckit::CodeLocation& loc) {
    throw eckit::NotImplemented(loc);
}

void throw_NotImplemented(const std::string& msg, const eckit::CodeLocation& loc) {
    throw eckit::NotImplemented(msg, loc);
}

void throw_AssertionFailed(const std::string& msg) {
    throw eckit::AssertionFailed(msg);
}

void throw_AssertionFailed(const std::string& msg, const eckit::CodeLocation& loc) {
    throw eckit::AssertionFailed(msg, loc);
}

void throw_AssertionFailed(const std::string& code, const std::string& msg, const eckit::CodeLocation& loc) {
    std::ostringstream ss;
    ss << " [[ " << code << " ]]\n" << msg;
    throw eckit::AssertionFailed(ss.str(), loc);
}

void throw_Exception(const std::string& msg) {
    throw eckit::Exception(msg);
}

void throw_Exception(const std::string& msg, const eckit::CodeLocation& loc) {
    throw eckit::Exception(msg, loc);
}

void throw_CantOpenFile(const std::string& file) {
    throw eckit::CantOpenFile(file);
}

void throw_CantOpenFile(const std::string& file, const eckit::CodeLocation& loc) {
    throw eckit::CantOpenFile(file, loc);
}

void throw_OutOfRange(const std::string& varname, idx_t index, idx_t size) {
    std::ostringstream ss;
    ss << "OutOfRange: Tried to access " << varname << " index " << index << " but maximum allowed index is "
       << size - 1;
    throw eckit::Exception(ss.str());
}

void throw_OutOfRange(const std::string& varname, idx_t index, idx_t size, const eckit::CodeLocation& loc) {
    std::ostringstream ss;
    ss << "OutOfRange: Tried to access " << varname << " index " << index << " but maximum allowed index is "
       << size - 1;
    throw eckit::Exception(ss.str(), loc);
}


}  // namespace atlas
