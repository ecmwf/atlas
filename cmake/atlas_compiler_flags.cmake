# (C) Copyright 1996-2014 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.


if( NOT CMAKE_BUILD_TYPE MATCHES "Debug" )
  add_definitions( -DNDEBUG )
endif( )

#######################################################################################
# Fortran
#######################################################################################

if( CMAKE_Fortran_COMPILER_ID MATCHES "GNU" )
  include( compiler_flags_GNU_Fortran )
elseif( CMAKE_Fortran_COMPILER_ID MATCHES "Intel" )
  include( compiler_flags_Intel_Fortran )
elseif( CMAKE_Fortran_COMPILER_ID MATCHES "XL" )
  include( compiler_flags_XL_Fortran )
elseif( CMAKE_Fortran_COMPILER_ID MATCHES "Cray" )
  include( compiler_flags_Cray_Fortran )
else()
  message( STATUS "Fortran compiler with ID ${CMAKE_Fortran_COMPILER_ID} will be used with CMake default options")
endif()

#######################################################################################
# C
#######################################################################################

# todo

#######################################################################################
# C++
#######################################################################################

if( CMAKE_CXX_COMPILER_ID MATCHES "GNU" )
  include( compiler_flags_GNU_CXX )
elseif( CMAKE_CXX_COMPILER_ID MATCHES "Intel" )
  include( compiler_flags_Intel_CXX )
elseif( CMAKE_CXX_COMPILER_ID MATCHES "XL" )
  include( compiler_flags_XL_CXX )
elseif( CMAKE_CXX_COMPILER_ID MATCHES "Cray" )
  include( compiler_flags_Cray_CXX )
elseif( CMAKE_CXX_COMPILER_ID MATCHES "Clang" )
  include( compiler_flags_Clang_CXX )
else()
  message( STATUS "C++ compiler with ID ${CMAKE_CXX_COMPILER_ID} will be used with CMake default options")
endif()
