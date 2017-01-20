#ifndef atlas_runtime_Log_h
#define atlas_runtime_Log_h

#include "atlas/internals/atlas_defines.h"

#ifdef ATLAS_HAVE_FORTRAN
#include "fckit/Log.h"
namespace atlas { 
typedef fckit::Log Log;
}

#else

#include "eckit/log/Log.h"
namespace atlas {
class Log : public eckit::Log {
  enum Style {
    SIMPLE=0,PREFIX=1,TIMESTAMP=2
  };
  static void addFortranUnit(int unit, Style=PREFIX, const char* prefix="") { /*NOTIMP*/ }
  static void setFortranUnit(int unit, Style=PREFIX, const char* prefix="") { /*NOTIMP*/ }

  // Fortran unit numbers
  static int output_unit() { return 6; }
  static int error_unit()  { return 0; }
};
}
#endif

#define ATLAS_DEBUG(WHAT)    do{ atlas::Log::info() << "DEBUG(" << WHAT << ") @ " << Here() << std::endl; } while(0)
#define ATLAS_DEBUG_HERE()   do{ atlas::Log::info() << "DEBUG() @ " << Here() << std::endl; } while(0)
#define ATLAS_DEBUG_VAR(VAR) do{ atlas::Log::info() << "DEBUG( " << #VAR << " : " << VAR << " ) @ " << Here() << std::endl; } while(0)

#endif
