#ifndef atlas_runtime_Log_h
#define atlas_runtime_Log_h

#include "eckit/log/Log.h"

namespace atlas {
  typedef eckit::Log Log;
} // namespace atlas

#define ATLAS_DEBUG(WHAT)    do{ atlas::Log::info() << "DEBUG(" << WHAT << ") @ " << Here() << std::endl; } while(0)
#define ATLAS_DEBUG_HERE()   do{ atlas::Log::info() << "DEBUG() @ " << Here() << std::endl; } while(0)
#define ATLAS_DEBUG_VAR(VAR) do{ atlas::Log::info() << "DEBUG( " << #VAR << " : " << VAR << " ) @ " << Here() << std::endl; } while(0)

#endif
