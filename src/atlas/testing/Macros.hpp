/*
 * (C) Copyright 1996-2014 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

/// @file Macros.hpp
/// @author Willem Deconinck
/// This file implements Macros to be used inside the UnitTest framework

#include <stdexcept>
#include <sstream>

// ------------------------------------------------------------------------------------------------

/*
Basic Macros and Defines to be used as building blocks
*/  

#define ATLAS_CHECK_SUCCESS 0
#define ATLAS_CHECK_FAILED  1

#define __ATLAS_CHECK_FATAL(CODE) { int line=__LINE__; const char* file=__FILE__; const char* func=__FUNCTION__; \
  int ierr = CODE; \
  if (ierr != ATLAS_CHECK_SUCCESS) { std::stringstream ss; \
    ss << "FATAL CHECK FAILED in function '"<<func<<"' ("<<file<<":"<<line<<") : " << std::string(#CODE) << " " << error_msg_; \
    throw std::runtime_error( ss.str() ); \
  } }
  
#define __ATLAS_CHECK(CODE) { int line=__LINE__; const char* file=__FILE__; const char* func=__FUNCTION__; \
  int ierr = CODE; \
  if (ierr != ATLAS_CHECK_SUCCESS) { std::stringstream ss; \
    ss << "CHECK FAILED in function '"<<func<<"' ("<<file<<":"<<line<<") : " << std::string(#CODE) << " " << error_msg_; \
    std::cout << ss.str() << std::endl; \
  } }

// ------------------------------------------------------------------------------------------------

/*
Following macro's are only to be used inside a UnitTest member function
*/

#define ATLAS_CHECK( EXPR )            __ATLAS_CHECK( assert(  EXPR  ) )
#define ATLAS_CHECK_FATAL( EXPR )      __ATLAS_CHECK_FATAL( assert(  EXPR  ) )
#define ATLAS_CHECK_EQUAL(v1,v2)       __ATLAS_CHECK( assert_equal(v1,v2) )
#define ATLAS_CHECK_EQUAL_FATAL(v1,v2) __ATLAS_CHECK_FATAL( assert_equal(v1,v2) )
#define ATLAS_CHECK_NOTHROW( CODE ) \
  try { CODE; } \
  catch(...) { \
    ss << "CHECK_NOTHROW FAILED in function '"<<func<<"' ("<<file<<":"<<line<<") : " << std::string(#CODE); \
  } 
#define ATLAS_CHECK_NOTHROW_FATAL( CODE ) \
  try { CODE; } \
  catch(...) { \
    ss << "CHECK_NOTHROW FAILED in function '"<<func<<"' ("<<file<<":"<<line<<") : " << std::string(#CODE); \
    throw; \
  }
