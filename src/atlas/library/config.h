/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#pragma once

#include <cstddef>

#include "eckit/eckit.h"

#include "atlas/atlas_ecbuild_config.h"
#include "atlas/library/defines.h"

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#define DOXYGEN_HIDE(X) X
#endif

#define ATLAS_HAVE_TRACE 1
#define ATLAS_HAVE_TRACE_BARRIERS 1

namespace atlas {

/// @typedef gidx_t
/// Integer type for global indices
#if ATLAS_BITS_GLOBAL == 32
using gidx_t = int;
#else
using gidx_t = long;
#endif

/// @typedef idx_t
/// Integer type for indices in connectivity tables
#if (ATLAS_BITS_LOCAL == 32)
using idx_t = int;
#else
using idx_t  = long;
#endif

/// @typedef uidx_t
/// Integer type for unique indices
typedef gidx_t uidx_t;

#define ATLAS_ECKIT_VERSION_AT_LEAST(x, y, z) (ATLAS_ECKIT_VERSION_INT >= x * 10000 + y * 100 + z)

#if ATLAS_ECKIT_VERSION_AT_LEAST(1, 19, 0)
#define ATLAS_ECKIT_HAVE_ECKIT_585 1
#else
#define ATLAS_ECKIT_HAVE_ECKIT_585 0
#endif

}  // namespace atlas
