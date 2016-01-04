# (C) Copyright 1996-2015 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

if( NOT CMAKE_TOOLCHAIN_FILE )
  
  ecbuild_debug("No toolchain file present, attempt with built-in compiler flags")

  #######################################################################################
  # Fortran
  #######################################################################################
  
  include( compiler_flags_${CMAKE_Fortran_COMPILER_ID}_Fortran OPTIONAL 
           RESULT_VARIABLE INCLUDED_Fortran_FLAGS )
  if( NOT INCLUDED_Fortran_FLAGS )
      ecbuild_info("Fortran compiler with ID ${CMAKE_Fortran_COMPILER_ID} will be used with CMake default options")
  endif()

  #######################################################################################
  # C
  #######################################################################################

  # todo

  #######################################################################################
  # C++
  #######################################################################################

  include( compiler_flags_${CMAKE_CXX_COMPILER_ID}_CXX OPTIONAL 
           RESULT_VARIABLE INCLUDED_CXX_FLAGS )
  if( NOT INCLUDED_CXX_FLAGS )
      ecbuild_info("C++ compiler with ID ${CMAKE_CXX_COMPILER_ID} will be used with CMake default options")
  endif()

else()
  
  ecbuild_debug("Toolchain file present, all compiler flags should already be present")
  
endif()
