/*
 * (C) Copyright 1996-2014 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef atlas_Debug_hpp
#define atlas_Debug_hpp

#include <sstream>
#include <unistd.h>

#include "atlas/mpl/MPL.hpp"
#include "atlas/atlas_config.h"


/// DEBUG MACRO
#define DEBUG_0()            std::cerr << "["<< MPL::rank() << "] DEBUG() @ " << Here() << std::endl;
#define DEBUG_1(WHAT)        std::cerr << "["<< MPL::rank() << "] DEBUG( " << WHAT << " ) @ " << Here() << std::endl;
#define DEBUG_2(WHAT,RANK)   if(MPL::rank() == RANK) { DEBUG_1(WHAT) }
#define DEBUG_X(x,A,B,FUNC, ...)  FUNC
#define DEBUG(...)  DEBUG_X(,##__VA_ARGS__,\
                        DEBUG_2(__VA_ARGS__),\
                        DEBUG_1(__VA_ARGS__),\
                        DEBUG_0(__VA_ARGS__))

/// DEBUG_SYNC MACRO
#define DEBUG_SYNC(...) \
  MPI_Barrier(MPI_COMM_WORLD);\
  DEBUG_X(,##__VA_ARGS__,\
     DEBUG_2(__VA_ARGS__),\
     DEBUG_1(__VA_ARGS__),\
     DEBUG_0(__VA_ARGS__))\
  MPI_Barrier(MPI_COMM_WORLD); usleep(1000); /*microseconds*/

/// DEBUG_VAR MACRO
#ifdef DEBUG_VAR
  #undef DEBUG_VAR
#endif
#define DEBUG_VAR_1(VAR) \
  std::cerr << "["<< MPL::rank() << "] DEBUG( " << #VAR << " : " << VAR << " ) @ " << Here() << std::endl;
#define DEBUG_VAR_2(VAR,RANK) if(MPL::rank() == RANK) { DEBUG_VAR_1(VAR) }
#define DEBUG_VAR_X(x,A,B,FUNC, ...)  FUNC
#define DEBUG_VAR(...)  DEBUG_VAR_X(,##__VA_ARGS__,\
                        DEBUG_VAR_2(__VA_ARGS__),\
                        DEBUG_VAR_1(__VA_ARGS__))

/// DEBUG_VAR_SYNC MACRO
#define DEBUG_VAR_SYNC(...) \
  MPI_Barrier(MPI_COMM_WORLD);\
  DEBUG_VAR_X(,##__VA_ARGS__,\
     DEBUG_VAR_2(__VA_ARGS__),\
     DEBUG_VAR_1(__VA_ARGS__))\
  MPI_Barrier(MPI_COMM_WORLD); usleep(1000); /*microseconds*/

#endif
