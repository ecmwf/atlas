#ifndef atlas_atlas_h
#define atlas_atlas_h

#include <iosfwd>

#include "atlas/internals/atlas_config.h"

// Forward declarations
namespace eckit
{
  class Parametrisation;
}

namespace atlas {

// This is only to be used from Fortran or unit-tests
void atlas_init(int argc, char **argv);

void atlas_init();
void atlas_init(const eckit::Parametrisation&);
void atlas_finalize();


void atlas_info(std::ostream&);

//----------------------------------------------------------------------------------------------------------------------

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
  const char* atlas__rundir ();
  const char* atlas__workdir ();
}

//----------------------------------------------------------------------------------------------------------------------

} // namespace atlas

#endif
