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

#include "atlas/library/config.h"
#include "eckit/log/CodeLocation.h"

namespace atlas {

[[noreturn]] void throw_NotImplemented( const eckit::CodeLocation& );
[[noreturn]] void throw_NotImplemented( const std::string&, const eckit::CodeLocation& );

[[noreturn]] void throw_AssertionFailed( const std::string& );
[[noreturn]] void throw_AssertionFailed( const std::string&, const eckit::CodeLocation& );

[[noreturn]] void throw_Exception( const std::string& );
[[noreturn]] void throw_Exception( const std::string&, const eckit::CodeLocation& );

[[noreturn]] void throw_SeriousBug( const std::string& );
[[noreturn]] void throw_SeriousBug( const std::string&, const eckit::CodeLocation& );

namespace detail {
inline void Assert( bool success, const char* msg, const char* file, int line, const char* func ) {
    if ( not success ) { throw_AssertionFailed( msg, eckit::CodeLocation( file, line, func ) ); }
}
}  // namespace detail
#define ATLAS_ASSERT( a ) ::atlas::detail::Assert( a, #a, __FILE__, __LINE__, __func__ )
#define ATLAS_NOTIMPLEMENTED ::atlas::throw_NotImplemented( Here() )

}  // namespace atlas
