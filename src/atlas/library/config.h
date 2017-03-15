#pragma once

#include "atlas/atlas_ecbuild_config.h"
#include "atlas/library/defines.h"

namespace atlas {

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
