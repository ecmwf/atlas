# (C) Copyright 1996-2014 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )
include_directories( ${CMAKE_CURRENT_BINARY_DIR} )

### mpistubs

debug_var( HAVE_MPI )

if( NOT HAVE_MPI )

	include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/mpistubs )
	list( APPEND TPL mpistubs )

	add_subdirectory( mpistubs )

endif()

### sources

add_subdirectory( atlas )

add_subdirectory( apps )

add_subdirectory( tests )

add_subdirectory( sandbox )
