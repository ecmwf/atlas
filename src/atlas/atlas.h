#ifndef atlas_atlas_h
#define atlas_atlas_h

#include "atlas/atlas_config.h"

namespace atlas {

void atlas_init(int argc=0, char **argv=0);
void atlas_finalize();

// ------------------------------------------------------------------
// C wrapper interfaces to C++ routines
extern "C"
{
  void atlas__atlas_init_noargs();
  void atlas__atlas_init (int argc, char** argv);
  void atlas__atlas_finalize ();
  const char* atlas__eckit_version();
  const char* atlas__eckit_git_sha1();
  const char* atlas__eckit_git_sha1_abbrev (int length);
  const char* atlas__atlas_version();
  const char* atlas__atlas_git_sha1();
  const char* atlas__atlas_git_sha1_abbrev (int length);
  const char* atlas__run_name ();
  const char* atlas__display_name ();
}
// ------------------------------------------------------------------


} // namespace atlas

#endif
