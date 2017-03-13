#pragma once

#include "atlas/atlas_ecbuild_config.h"
#include "atlas/library/defines.h"

namespace atlas {
  
  const char* version();

  // VVMMPP ( 10000*ATLAS_MAJOR_VERSION + 100*ATLAS_MINOR_VERSION + ATLAS_PATCH_VERSION )
  int version_int();

  const char* git_sha1(unsigned int chars=7);

/// @typedef gidx_t
/// Integer type for global indices
#if ATLAS_BITS_GLOBAL==32
  typedef int  gidx_t;
#else
  typedef long gidx_t;
#endif

  /// @typedef idx_t
  /// Integer type for indices in connectivity tables
  typedef int idx_t;

  /// @typedef uidx_t
  /// Integer type for unique indices
  typedef long uidx_t;

}
