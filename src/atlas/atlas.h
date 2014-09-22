#ifndef atlas_atlas_h
#define atlas_atlas_h

#include "atlas/atlas_config.h"

namespace atlas {

void atlas_init(int argc, char** argv);
void atlas_finalize();

// ------------------------------------------------------------------
// C wrapper interfaces to C++ routines
extern "C"
{
  void atlas__atlas_init (int argc, char** argv);
	void atlas__atlas_finalize ();
}
// ------------------------------------------------------------------


} // namespace atlas

#endif
