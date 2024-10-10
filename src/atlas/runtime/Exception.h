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
#include <sstream>

#include "eckit/log/CodeLocation.h"

#include "atlas/library/config.h"

namespace atlas {

#ifndef DOXYGEN_SHOULD_SKIP_THIS
[[noreturn]] void throw_NotImplemented(const eckit::CodeLocation&);
[[noreturn]] void throw_NotImplemented(const std::string&, const eckit::CodeLocation&);

[[noreturn]] void throw_AssertionFailed(const std::string&);
[[noreturn]] void throw_AssertionFailed(const std::string&, const eckit::CodeLocation&);
[[noreturn]] void throw_AssertionFailed(const std::string&, const std::string&, const eckit::CodeLocation&);

[[noreturn]] void throw_Exception(const std::string&);
[[noreturn]] void throw_Exception(const std::string&, const eckit::CodeLocation&);

[[noreturn]] void throw_CantOpenFile(const std::string&);
[[noreturn]] void throw_CantOpenFile(const std::string&, const eckit::CodeLocation&);

[[noreturn]] void throw_OutOfRange(const std::string& varname, idx_t index, idx_t size);
[[noreturn]] void throw_OutOfRange(const std::string& varname, idx_t index, idx_t size, const eckit::CodeLocation&);
[[noreturn]] void throw_OutOfRange(idx_t index, idx_t size);
[[noreturn]] void throw_OutOfRange(idx_t index, idx_t size, const eckit::CodeLocation&);

namespace detail {
inline void Assert(bool success, const char* code, const char* file, int line, const char* func) {
    if (not success) {
        throw_AssertionFailed(code, eckit::CodeLocation(file, line, func));
    }
}
inline void Assert(bool success, const char* code, const std::string& msg, const char* file, int line,
                   const char* func) {
    if (not success) {
        throw_AssertionFailed(code, msg, eckit::CodeLocation(file, line, func));
    }
}

}  // namespace detail
}  // namespace atlas

#include "atlas/library/detail/BlackMagic.h"

#define ATLAS_NOTIMPLEMENTED ::atlas::throw_NotImplemented(Here())


#define ATLAS_ASSERT_NOMSG(a) \
    static_cast<void>(0), (a) ? (void)0 : ::atlas::detail::Assert(bool(a), #a, __FILE__, __LINE__, __func__)
#define ATLAS_ASSERT_MSG(a, m) \
    static_cast<void>(0), (a) ? (void)0 : ::atlas::detail::Assert(bool(a), #a, m, __FILE__, __LINE__, __func__)

#define ATLAS_ASSERT(...) __ATLAS_SPLICE(__ATLAS_ASSERT_, __ATLAS_NARG(__VA_ARGS__))(__VA_ARGS__)
#define __ATLAS_ASSERT_1 ATLAS_ASSERT_NOMSG
#define __ATLAS_ASSERT_2 ATLAS_ASSERT_MSG

#define ATLAS_THROW_EXCEPTION(WHAT)                 \
    do {                                            \
        std::ostringstream ss;                      \
        ss << WHAT;                                 \
        ::atlas::throw_Exception(ss.str(), Here()); \
    } while (false)

#endif
