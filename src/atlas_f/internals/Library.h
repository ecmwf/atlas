#pragma once

// C wrapper interfaces to C++ routines
extern "C"
{
  void atlas__atlas_init_noargs();
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
