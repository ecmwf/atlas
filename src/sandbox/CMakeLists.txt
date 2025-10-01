# (C) Copyright 2013 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

if( ECKIT_INCLUDE_DIRS ) # eckit not yet ported to CMake3
    include_directories( ${ECKIT_INCLUDE_DIRS} )
endif()

add_subdirectory( apps )
if( HAVE_FORTRAN )
add_subdirectory( fortran_submodule )
add_subdirectory( fortran_object )
add_subdirectory( example_fortran )
add_subdirectory( fortran_acc_fields )
endif()
add_subdirectory( interpolation )
add_subdirectory( interpolation-fortran )
add_subdirectory( grid_distribution )
add_subdirectory( benchmark_ifs_setup )
add_subdirectory( benchmark_sorting )
add_subdirectory( benchmark_trans )
add_subdirectory( read_field_api_output )
