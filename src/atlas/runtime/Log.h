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

#include "atlas/library/config.h"

#if ATLAS_HAVE_FORTRAN
#include "fckit/Log.h"
namespace atlas {
namespace detail {
typedef fckit::Log LogBase;
}
}  // namespace atlas
#else
#include "eckit/log/Log.h"
namespace atlas {
namespace detail {
typedef eckit::Log LogBase;
}
}  // namespace atlas
#endif

namespace atlas {

class Log : public detail::LogBase {
public:
    using Channel = eckit::Channel;  // derives from std::ostream

    static Channel& info();
    static Channel& trace();
    static Channel& debug();

#if !ATLAS_HAVE_FORTRAN
    // Stubs for what fckit::Log provides
    enum Style
    {
        SIMPLE    = 0,
        PREFIX    = 1,
        TIMESTAMP = 2
    };
    static void addFortranUnit( int unit, Style = PREFIX, const char* prefix = "" ) { /*NOTIMP*/
    }
    static void setFortranUnit( int unit, Style = PREFIX, const char* prefix = "" ) { /*NOTIMP*/
    }

    // Fortran unit numbers
    static int output_unit() { return 6; }
    static int error_unit() { return 0; }
#endif
};

std::string backtrace();

// Non-flushing version of std::endl
inline std::ostream& newl( std::ostream& out ) {
    return out << '\n';
}

}  // namespace atlas

#include <sstream>
#include "atlas/util/detail/BlackMagic.h"
#include "eckit/log/CodeLocation.h"

namespace atlas {
namespace detail {
void debug_parallel_here( const eckit::CodeLocation& );
void debug_parallel_what( const eckit::CodeLocation&, const std::string& );
}  // namespace detail
}  // namespace atlas

#define ATLAS_DEBUG_HERE()                                           \
    do {                                                             \
        ::atlas::Log::info() << "DEBUG() @ " << Here() << std::endl; \
    } while ( 0 )
#define ATLAS_DEBUG_WHAT( WHAT )                                                   \
    do {                                                                           \
        ::atlas::Log::info() << "DEBUG(" << WHAT << ") @ " << Here() << std::endl; \
    } while ( 0 )
#define ATLAS_DEBUG_VAR( VAR )                                                                       \
    do {                                                                                             \
        ::atlas::Log::info() << "DEBUG( " << #VAR << " : " << VAR << " ) @ " << Here() << std::endl; \
    } while ( 0 )

#define ATLAS_DEBUG( ... ) __ATLAS_SPLICE( __ATLAS_DEBUG_, __ATLAS_NARG( __VA_ARGS__ ) )( __VA_ARGS__ )
#define __ATLAS_DEBUG_0 ATLAS_DEBUG_HERE
#define __ATLAS_DEBUG_1 ATLAS_DEBUG_WHAT

#define ATLAS_DEBUG_BACKTRACE()                                                                                \
    do {                                                                                                       \
        ::atlas::Log::info() << "DEBUG() @ " << Here() << "Backtrace:\n" << ::atlas::backtrace() << std::endl; \
    } while ( 0 )

#define ATLAS_DEBUG_PARALLEL_HERE()                     \
    do {                                                \
        ::atlas::detail::debug_parallel_here( Here() ); \
    } while ( 0 )
#define ATLAS_DEBUG_PARALLEL_WHAT( WHAT )                        \
    do {                                                         \
        std::stringstream w;                                     \
        w << WHAT;                                               \
        ::atlas::detail::debug_parallel_what( Here(), w.str() ); \
    } while ( 0 )

#define ATLAS_DEBUG_PARALLEL( ... )                                        \
    __ATLAS_SPLICE( __ATLAS_DEBUG_PARALLEL_, __ATLAS_NARG( __VA_ARGS__ ) ) \
    ( __VA_ARGS__ )
#define __ATLAS_DEBUG_PARALLEL_0 ATLAS_DEBUG_PARALLEL_HERE
#define __ATLAS_DEBUG_PARALLEL_1 ATLAS_DEBUG_PARALLEL_WHAT
