# (C) Copyright 1996-2016 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

# - Try to find gridtools_storage
# Once done this will define
#  GRIDTOOLS_STORAGE_FOUND - True if gridtools_storage found
#  GRIDTOOLS_STORAGE_INCLUDE_DIRS - The gridtools_storage include directories
#  GRIDTOOLS_STORAGE_LIBRARIES - The libraries needed to use gridtools_storage
#  GRIDTOOLS_STORAGE_DEFINITIONS - Compiler switches required for using gridtools_storage

if( NOT GRIDTOOLS_STORAGE_FOUND )

  find_path( GRIDTOOLS_STORAGE_INCLUDE_DIR 
             NAMES storage/storage-facility.hpp
             PATHS 
                ${CMAKE_INSTALL_PREFIX}
                "${GRIDTOOLS_STORAGE_PATH}"
                ENV GRIDTOOLS_STORAGE_PATH 
            PATH_SUFFIXES include
  )
  
  ecbuild_debug( "Searching for Boost version >= 1.58, required for gridtools_storage..." )
  if( GRIDTOOLS_STORAGE_INCLUDE_DIR AND NOT Boost_FOUND )
    find_package(Boost 1.58.0 )
  endif()
  if( Boost_FOUND )
    ecbuild_debug( "Boost, required for gridtools_storage found at ${Boost_INCLUDE_DIRS}" )
  else()
    ecbuild_debug( "Boost, required for gridtools_storage not found" )
  endif()
  
  include(FindPackageHandleStandardArgs)

  # handle the QUIETLY and REQUIRED arguments and set GRIDTOOLS_STORAGE_FOUND to TRUE
  find_package_handle_standard_args( gridtools_storage DEFAULT_MSG
                                     GRIDTOOLS_STORAGE_INCLUDE_DIR Boost_INCLUDE_DIRS )

  mark_as_advanced( GRIDTOOLS_STORAGE_INCLUDE_DIRS GRIDTOOLS_STORAGE_LIBRARIES )

  set( gridtools_storage_FOUND ${GRIDTOOLS_STORAGE_FOUND} )
  set( GRIDTOOLS_STORAGE_INCLUDE_DIRS 
    ${GRIDTOOLS_STORAGE_INCLUDE_DIR}
    ${Boost_INCLUDE_DIRS}  
  )

  # message("DEBUG Findgridtools_storage : GRIDTOOLS_STORAGE_FOUND         [${GRIDTOOLS_STORAGE_FOUND}]")
  # message("DEBUG Findgridtools_storage : gridtools_storage_FOUND         [${gridtools_storage_FOUND}]")
  # message("DEBUG Findgridtools_storage : GRIDTOOLS_STORAGE_INCLUDE_DIRS  [${GRIDTOOLS_STORAGE_INCLUDE_DIRS}]")

endif()
