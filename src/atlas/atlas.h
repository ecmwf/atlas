#ifndef atlas_atlas_h
#define atlas_atlas_h

#include "atlas_config.h"
#include "atlas_defines.h"
#include "atlas_version.h"

namespace atlas {

void atlas_init(int argc, char** argv);

// ------------------------------------------------------------------
// C wrapper interfaces to C++ routines
extern "C"
{
  void atlas__atlas_init (int argc, char** argv);
}
// ------------------------------------------------------------------


} // namespace atlas

#endif
