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

#include <iosfwd>
#include <string>
#include "eckit/log/CodeLocation.h"

namespace atlas {

class CodeLocation {
public:
    CodeLocation(const CodeLocation& loc): CodeLocation(loc.file(), loc.line(), loc.func(), loc.stored_) {}
    CodeLocation(const eckit::CodeLocation& loc): loc_(loc), stored_(false) {}
    CodeLocation(const char* file, int line, const char* function, bool store = false): stored_(store) {
        if (stored_) {
            if (file) {
                file_str_ = std::string(file);
                file_     = file_str_.c_str();
            }
            if (function) {
                function_str_ = std::string(function);
                function_     = function_str_.c_str();
            }
            loc_ = eckit::CodeLocation(file_, line, function_);
        }
        else {
            loc_ = eckit::CodeLocation(file, line, function);
        }
    }
    operator const eckit::CodeLocation&() const { return loc_; }

    /// conversion to bool for checking if location was set
    operator bool() const { return loc_; }

    std::string asString() const { return loc_.asString(); }
    /// accessor to line
    int line() const { return loc_.line(); }
    /// accessor to file
    const char* file() const { return loc_.file(); }
    /// accessor to function
    const char* func() const { return loc_.func(); }
    friend std::ostream& operator<<(std::ostream& s, const CodeLocation& loc);

private:
    eckit::CodeLocation loc_;
    bool stored_          = false;
    const char* file_     = nullptr;
    const char* function_ = nullptr;
    std::string file_str_;
    std::string function_str_;
};

}  // namespace atlas
