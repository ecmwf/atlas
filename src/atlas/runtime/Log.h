#pragma once

#include "atlas/library/config.h"

#include "atlas/library/Library.h"

#ifdef ATLAS_HAVE_FORTRAN
#include "fckit/Log.h"
namespace atlas { namespace detail { typedef fckit::Log LogBase; } }
#else
#include "eckit/log/Log.h"
namespace atlas { namespace detail { typedef eckit::Log LogBase; } }
#endif

namespace atlas {

class Log : public detail::LogBase {

public:

    typedef eckit::Channel Channel;

    static std::ostream& info()  { return atlas::Library::instance().infoChannel();  }
    static std::ostream& trace() { return atlas::Library::instance().traceChannel(); }
    static std::ostream& debug() { return atlas::Library::instance().debugChannel(); }

    static std::ostream& debug_parallel();

#ifndef ATLAS_HAVE_FORTRAN
    // Stubs for what fckit::Log provides
    enum Style {
      SIMPLE=0,PREFIX=1,TIMESTAMP=2
    };
    static void addFortranUnit(int unit, Style=PREFIX, const char* prefix="") { /*NOTIMP*/ }
    static void setFortranUnit(int unit, Style=PREFIX, const char* prefix="") { /*NOTIMP*/ }

    // Fortran unit numbers
    static int output_unit() { return 6; }
    static int error_unit()  { return 0; }
#endif
};

std::string backtrace();

namespace detail {
  void print_parallel_here(std::ostream&,const eckit::CodeLocation&);
}

} // namespace atlas

#include "atlas/util/detail/BlackMagic.h"

#define ATLAS_DEBUG_HERE()     do{ ::atlas::Log::info() << "DEBUG() @ " << Here() << std::endl; } while(0)
#define ATLAS_DEBUG_WHAT(WHAT) do{ ::atlas::Log::info() << "DEBUG(" << WHAT << ") @ " << Here() << std::endl; } while(0)
#define ATLAS_DEBUG_VAR(VAR)   do{ ::atlas::Log::info() << "DEBUG( " << #VAR << " : " << VAR << " ) @ " << Here() << std::endl; } while(0)

#define ATLAS_DEBUG(...) __ATLAS_SPLICE( __ATLAS_DEBUG_, __ATLAS_NARGS(__VA_ARGS__) ) (__VA_ARGS__)
#define __ATLAS_DEBUG_0 ATLAS_DEBUG_HERE
#define __ATLAS_DEBUG_1 ATLAS_DEBUG_WHAT

#define ATLAS_DEBUG_BACKTRACE() do{ ::atlas::Log::info() << "DEBUG() @ " << Here() << "Backtrace:\n" << backtrace() << std::endl; } while (0)

#define ATLAS_DEBUG_PARALLEL_HERE() do{ ::atlas::print_parallel_here(std::cout,Here()); } while (0)
