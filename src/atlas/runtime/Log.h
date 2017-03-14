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

// Use Log::debug<Atlas>() for debugging that can be switched
// on or off by the boolean environment variable ATLAS_DEBUG

class Log : public detail::LogBase {

public:

    typedef eckit::Channel Channel;

#ifndef ATLAS_HAVE_FORTRAN
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

} // namespace atlas
