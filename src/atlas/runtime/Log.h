#pragma once

#include "atlas/library/config.h"

#include "eckit/log/Log.h"

namespace atlas {

class Log : public eckit::Log {
public:
    using Channel = eckit::Channel;  // derives from std::ostream

    static Channel& info();
    static Channel& warning();
    static Channel& trace();
    static Channel& debug();

    // Same as what fckit::Log provides
    enum Style
    {
        SIMPLE    = 0,
        PREFIX    = 1,
        TIMESTAMP = 2
    };
    static void addFortranUnit(int unit, Style = PREFIX, const char* prefix = "");
    static void setFortranUnit(int unit, Style = PREFIX, const char* prefix = "");

    // Fortran unit numbers
    static int output_unit() { return 6; }
    static int error_unit() { return 0; }
};

std::string backtrace();

}  // namespace atlas

#include <sstream>
#include "atlas/library/detail/BlackMagic.h"
#include "eckit/log/CodeLocation.h"

namespace atlas {
namespace detail {
void debug_parallel_here(const eckit::CodeLocation&);
void debug_parallel_what(const eckit::CodeLocation&, const std::string&);
}  // namespace detail
}  // namespace atlas

#define ATLAS_DEBUG_HERE()                                           \
    do {                                                             \
        ::atlas::Log::info() << "DEBUG() @ " << Here() << std::endl; \
    } while (0)
#define ATLAS_DEBUG_WHAT(WHAT)                                                     \
    do {                                                                           \
        ::atlas::Log::info() << "DEBUG(" << WHAT << ") @ " << Here() << std::endl; \
    } while (0)
#define ATLAS_DEBUG_VAR(VAR)                                                                         \
    do {                                                                                             \
        ::atlas::Log::info() << "DEBUG( " << #VAR << " : " << VAR << " ) @ " << Here() << std::endl; \
    } while (0)

#define ATLAS_DEBUG(...) __ATLAS_SPLICE(__ATLAS_DEBUG_, __ATLAS_NARG(__VA_ARGS__))(__VA_ARGS__)
#define __ATLAS_DEBUG_0 ATLAS_DEBUG_HERE
#define __ATLAS_DEBUG_1 ATLAS_DEBUG_WHAT

#define ATLAS_DEBUG_BACKTRACE()                                                                                \
    do {                                                                                                       \
        ::atlas::Log::info() << "DEBUG() @ " << Here() << "Backtrace:\n" << ::atlas::backtrace() << std::endl; \
    } while (0)

#define ATLAS_DEBUG_PARALLEL_HERE()                   \
    do {                                              \
        ::atlas::detail::debug_parallel_here(Here()); \
    } while (0)
#define ATLAS_DEBUG_PARALLEL_WHAT(WHAT)                        \
    do {                                                       \
        std::stringstream w;                                   \
        w << WHAT;                                             \
        ::atlas::detail::debug_parallel_what(Here(), w.str()); \
    } while (0)

#define ATLAS_DEBUG_PARALLEL(...)                                      \
    __ATLAS_SPLICE(__ATLAS_DEBUG_PARALLEL_, __ATLAS_NARG(__VA_ARGS__)) \
    (__VA_ARGS__)
#define __ATLAS_DEBUG_PARALLEL_0 ATLAS_DEBUG_PARALLEL_HERE
#define __ATLAS_DEBUG_PARALLEL_1 ATLAS_DEBUG_PARALLEL_WHAT
