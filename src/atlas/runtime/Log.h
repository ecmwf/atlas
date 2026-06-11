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

namespace eckit {
namespace mpi {
class Comm;
}  // namespace mpi
}  // namespace eckit

#include <sstream>
#include <string_view>
#include "atlas/library/detail/BlackMagic.h"
#include "eckit/log/CodeLocation.h"

namespace atlas {
namespace detail {
void debug_sync(const eckit::CodeLocation&, const eckit::mpi::Comm&);
void debug_sync(const eckit::CodeLocation&);
void debug_sync(const eckit::CodeLocation&, const eckit::mpi::Comm&, std::string_view);
void debug_sync(const eckit::CodeLocation&, std::string_view);
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

#define ATLAS_DEBUG_SYNC_HERE()                \
    do {                                           \
        ::atlas::detail::debug_sync(Here());   \
    } while (0)
#define ATLAS_DEBUG_SYNC_HERE_1(ARG1)              \
    do {                                           \
        ::atlas::detail::debug_sync(Here(), ARG1); \
    } while (0)
#define ATLAS_DEBUG_SYNC_HERE_2(ARG1, ARG2)              \
    do {                                                \
        ::atlas::detail::debug_sync(Here(), ARG1, ARG2); \
    } while (0)

#define ATLAS_DEBUG_SYNC(...)                                      \
    __ATLAS_SPLICE(__ATLAS_DEBUG_SYNC_, __ATLAS_NARG(__VA_ARGS__)) \
    (__VA_ARGS__)
#define __ATLAS_DEBUG_SYNC_0 ATLAS_DEBUG_SYNC_HERE
#define __ATLAS_DEBUG_SYNC_1 ATLAS_DEBUG_SYNC_HERE_1
#define __ATLAS_DEBUG_SYNC_2 ATLAS_DEBUG_SYNC_HERE_2

#define __ATLAS_DEBUG_SYNC_VAR_1(VAR)                                                                 \
    do {                                                                                              \
        std::ostringstream oss_debug_var;                                                             \
        oss_debug_var << #VAR << " : " << VAR;                                                        \
        ::atlas::detail::debug_sync(Here(), oss_debug_var.str());                                     \
    } while (0)
#define __ATLAS_DEBUG_SYNC_VAR_2(COMM, VAR)                                                          \
    do {                                                                                              \
        std::ostringstream oss_debug_var;                                                             \
        oss_debug_var << #VAR << " : " << VAR;                                                        \
        ::atlas::detail::debug_sync(Here(), COMM, oss_debug_var.str());                               \
    } while (0)

#define ATLAS_DEBUG_SYNC_VAR(...)                                        \
    __ATLAS_SPLICE(__ATLAS_DEBUG_SYNC_VAR_, __ATLAS_NARG(__VA_ARGS__)) \
    (__VA_ARGS__)
