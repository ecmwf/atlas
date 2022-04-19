# (C) Copyright 2011- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

# - Try to find the proj library
# Once done this will define
#
# PROJ_FOUND - system has proj
# PROJ_INCLUDE_DIRS - the proj include directory
# PROJ_LIBRARIES - Link these to use proj
# PROJ_VERSION - semantic version of proj
#

### Set search paths from environment
if ( NOT PROJ_PATH AND PROJ_ROOT )
    set( PROJ_PATH ${PROJ_ROOT} )
endif()
if ( NOT PROJ_PATH AND NOT "$ENV{PROJ_ROOT}" STREQUAL "" )
    set( PROJ_PATH "$ENV{PROJ_ROOT}" )
endif()
if ( NOT PROJ_PATH AND NOT "$ENV{PROJ_PATH}" STREQUAL "" )
    set( PROJ_PATH "$ENV{PROJ_PATH}" )
endif()
if ( NOT PROJ_PATH AND NOT "$ENV{PROJ_DIR}" STREQUAL "" )
    set( PROJ_PATH "$ENV{PROJ_DIR}" )
endif()


### If search paths given, use it, otherwise, use pkg-config
if( PROJ_PATH )
    find_path(PROJ_INCLUDE_DIR NAMES proj.h PATHS ${PROJ_PATH} ${PROJ_PATH}/include PATH_SUFFIXES PROJ NO_DEFAULT_PATH )
    find_library(PROJ_LIBRARY  NAMES proj   PATHS ${PROJ_PATH} ${PROJ_PATH}/lib     PATH_SUFFIXES PROJ NO_DEFAULT_PATH )
else()
    find_package(PkgConfig)
    if(PKG_CONFIG_FOUND)
        if(PROJ_FIND_VERSION)
            pkg_check_modules(PKPROJ ${_pkgconfig_REQUIRED} QUIET PROJ>=${PROJ_FIND_VERSION})
        else()
            pkg_check_modules(PKPROJ ${_pkgconfig_REQUIRED} QUIET PROJ)
        endif()

        if( PKPROJ_FOUND )
            find_path(PROJ_INCLUDE_DIR proj.h HINTS ${PKPROJ_INCLUDEDIR} ${PKPROJ_INCLUDE_DIRS} PATH_SUFFIXES PROJ NO_DEFAULT_PATH )
            find_library(PROJ_LIBRARY  proj   HINTS ${PKPROJ_LIBDIR}     ${PKPROJ_LIBRARY_DIRS} PATH_SUFFIXES PROJ NO_DEFAULT_PATH )
        endif()
    endif()
endif()

find_path(PROJ_INCLUDE_DIR NAMES proj.h PATHS PATH_SUFFIXES PROJ )
find_library( PROJ_LIBRARY NAMES proj   PATHS PATH_SUFFIXES PROJ )

### Detect version
set( PROJ_VERSION 0 )
file(READ ${PROJ_INCLUDE_DIR}/proj.h proj_version)
string(REGEX REPLACE "^.*PROJ_VERSION_MAJOR +([0-9]+).*$" "\\1" PROJ_VERSION_MAJOR "${proj_version}")
string(REGEX REPLACE "^.*PROJ_VERSION_MINOR +([0-9]+).*$" "\\1" PROJ_VERSION_MINOR "${proj_version}")
string(REGEX REPLACE "^.*PROJ_VERSION_PATCH +([0-9]+).*$" "\\1" PROJ_VERSION_PATCH "${proj_version}")
string(CONCAT PROJ_VERSION ${PROJ_VERSION_MAJOR} "." ${PROJ_VERSION_MINOR} "." ${PROJ_VERSION_PATCH})

### Handle the QUIETLY and REQUIRED arguments and set PROJ_FOUND
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(PROJ
                                  REQUIRED_VARS PROJ_LIBRARY PROJ_INCLUDE_DIR
                                  VERSION_VAR PROJ_VERSION)

set( PROJ_LIBRARIES    ${PROJ_LIBRARY} )
set( PROJ_INCLUDE_DIRS ${PROJ_INCLUDE_DIR} )

mark_as_advanced( PROJ_INCLUDE_DIR PROJ_LIBRARY )
