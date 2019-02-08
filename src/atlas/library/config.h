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
