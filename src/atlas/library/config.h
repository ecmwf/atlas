#pragma once

#include <cstddef>

#include "eckit/eckit.h"

#include "atlas/atlas_ecbuild_config.h"
#include "atlas/library/defines.h"

#define ATLAS_HAVE_TRACE 1

namespace atlas {

/// @typedef gidx_t
/// Integer type for global indices
#if ATLAS_BITS_GLOBAL == 32
typedef int gidx_t;
#else
typedef long gidx_t;
#endif

/// @typedef idx_t
/// Integer type for indices in connectivity tables
#if ( ATLAS_BITS_LOCAL == 32 )
typedef int idx_t;
#else
typedef long idx_t;
#endif

/// @typedef uidx_t
/// Integer type for unique indices
typedef gidx_t uidx_t;
}  // namespace atlas
