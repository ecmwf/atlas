/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef atlas_Debug_h
#define atlas_Debug_h

#include <sstream>
#include <unistd.h>
#include "atlas/internals/atlas_config.h"
#include "atlas/runtime/Log.h"
#include "atlas/parallel/mpi/mpi.h"

/// DEBUG MACRO
#define DEBUG_RANK (eckit::mpi::initialized() ? eckit::mpi::rank() : 0)
#define DEBUG_HERE()            atlas::Log::info() << "["<< DEBUG_RANK << "] DEBUG() @ " << Here() << std::endl;
#define DEBUG(WHAT)        atlas::Log::info() << "["<< DEBUG_RANK << "] DEBUG( " << WHAT << " ) @ " << Here() << std::endl;

/// DEBUG_VAR MACRO
#ifdef DEBUG_VAR
  #undef DEBUG_VAR
#endif
#define DEBUG_VAR(VAR) \
  atlas::Log::info() << "["<< DEBUG_RANK << "] DEBUG( " << #VAR << " : " << VAR << " ) @ " << Here() << std::endl;

#endif
