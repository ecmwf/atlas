# (C) Copyright 2025- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

# - Try to find the pocketfft library
# Once done this will define
#
# pocketfft_FOUND - pocketfft found
# pocketfft_INCLUDE_DIRS - the pocketfft include directories
# pocketfft_VERSION - semantic version of pocketfft
#

### Set search paths from environment
if ( NOT pocketfft_PATH AND pocketfft_ROOT )
    set( pocketfft_PATH ${pocketfft_ROOT} )
endif()
if ( NOT pocketfft_PATH AND NOT "$ENV{pocketfft_ROOT}" STREQUAL "" )
    set( pocketfft_PATH "$ENV{pocketfft_ROOT}" )
endif()
if ( NOT pocketfft_PATH AND NOT "$ENV{pocketfft_PATH}" STREQUAL "" )
    set( pocketfft_PATH "$ENV{pocketfft_PATH}" )
endif()

### If search paths given, use it
if( pocketfft_PATH )
    find_path(pocketfft_INCLUDE_DIR NAMES pocketfft_hdronly.h PATHS ${pocketfft_PATH} ${pocketfft_PATH}/include pocketfft NO_DEFAULT_PATH )
else()
    if( NOT pocketfft_INCLUDE_DIR )
        find_path(pocketfft_INCLUDE_DIR NAMES pocketfft_hdronly.h )
    endif()
endif()

### As of 2025 pocketfft does not yet have any version.

### Handle the QUIETLY and REQUIRED arguments and set pocketfft_FOUND
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(pocketfft
                                  REQUIRED_VARS pocketfft_INCLUDE_DIR)

set( POCKETFFT_INCLUDE_DIRS ${pocketfft_INCLUDE_DIR} )

mark_as_advanced( pocketfft_INCLUDE_DIR )
