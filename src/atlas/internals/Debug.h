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
#define DEBUG_RANK (eckit::mpi::Environment::instance().initialized() ? eckit::mpi::comm().rank() : 0)
#define DEBUG_0()            atlas::Log::info() << "["<< DEBUG_RANK << "] DEBUG() @ " << Here() << std::endl;
#define DEBUG_1(WHAT)        atlas::Log::info() << "["<< DEBUG_RANK << "] DEBUG( " << WHAT << " ) @ " << Here() << std::endl;
#define DEBUG_2(WHAT,RANK)   if(DEBUG_RANK == RANK) { DEBUG_1(WHAT) }
#define DEBUG_X(x,A,B,FUNC, ...)  FUNC
#define DEBUG(...)  do {DEBUG_X(,##__VA_ARGS__,\
                        DEBUG_2(__VA_ARGS__),\
                        DEBUG_1(__VA_ARGS__),\
                        DEBUG_0(__VA_ARGS__))} while(0)

/// DEBUG_SYNC MACRO
#define DEBUG_SYNC(...) do {\
  {eckit::mpi::comm().barrier();\
  int npid = eckit::mpi::comm().size();\
  for( int p=0; p<npid; ++p )\
  {\
    if( p==eckit::mpi::comm().rank() )\
    {\
      DEBUG_X(,##__VA_ARGS__,\
        DEBUG_2(__VA_ARGS__),\
        DEBUG_1(__VA_ARGS__),\
        DEBUG_0(__VA_ARGS__))\
    }\
    eckit::mpi::comm().barrier(); usleep(100); /*microseconds*/ \
  }}} while(0)

/// DEBUG_VAR MACRO
#ifdef DEBUG_VAR
  #undef DEBUG_VAR
#endif
#define DEBUG_VAR_1(VAR) \
  atlas::Log::info() << "["<< DEBUG_RANK << "] DEBUG( " << #VAR << " : " << VAR << " ) @ " << Here() << std::endl;
#define DEBUG_VAR_2(VAR,RANK) if(DEBUG_RANK == RANK) { DEBUG_VAR_1(VAR) }
#define DEBUG_VAR_X(x,A,B,FUNC, ...)  FUNC
#define DEBUG_VAR(...)  do {DEBUG_VAR_X(,##__VA_ARGS__,\
                            DEBUG_VAR_2(__VA_ARGS__),\
                            DEBUG_VAR_1(__VA_ARGS__))} while(0)

/// DEBUG_VAR_SYNC MACRO
#define DEBUG_VAR_SYNC(...) do {\
  eckit::mpi::comm().barrier();\
  DEBUG_VAR_X(,##__VA_ARGS__,\
     DEBUG_VAR_2(__VA_ARGS__),\
     DEBUG_VAR_1(__VA_ARGS__))\
  eckit::mpi::comm().barrier(); usleep(1000); /*microseconds*/\
  } while(0)

#endif
