# (C) Copyright 1996-2014 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

############################################################################################
# ATLAS
#
#		TODO:
#         * add support for triangulating spheres with CGAL
#

cmake_minimum_required( VERSION 2.8.4 FATAL_ERROR )

project( atlas Fortran CXX )

set( CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/../ecbuild/cmake")
set( ECBUILD_DEFAULT_BUILD_TYPE Release )

include( ecbuild_system )

ecbuild_requires_macro_version( 1.3 )

ecbuild_declare_project()

# fortran module support

set( CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules )
file( MAKE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY} )
include_directories( ${CMAKE_Fortran_MODULE_DIRECTORY} )

include( atlas_compiler_flags )

################################################################################################
# options & dependencies

ecbuild_find_python()

### OMP ...

ecbuild_add_option( FEATURE OMP
					DEFAULT ON
					DESCRIPTION "support for OpenMP threaded parallelism" )

if( HAVE_OMP )

  # !! TODO !! add test for OpenMP

  add_definitions( -DHAVE_OMP )

endif( )

### MPI ...

ecbuild_add_option( FEATURE MPI
					DEFAULT ON
					DESCRIPTION "support for MPI distributed parallelism" )

if( ENABLE_MPI )
   ecbuild_find_mpi()
   if( MPI_CXX_FOUND )
	  add_definitions( -DHAVE_MPI )
   else()
	  set( HAVE_MPI 0 )
   endif()
endif()

#get_filename_component(Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)
#if( NOT Fortran_COMPILER_NAME MATCHES "mp.*" AND NOT Fortran_COMPILER_NAME MATCHES "ftn" )
#  # message( FATAL_ERROR "MPI Compiler needed.\n"
#  #   "Please remove this build, export FC and CXX to a MPI compiler such as mpif90 or ftn on Cray and rerun." )
#  set( HAVE_MPI OFF )
#endif( )
#get_filename_component(CXX_COMPILER_NAME ${CMAKE_CXX_COMPILER} NAME)
#if( NOT CXX_COMPILER_NAME MATCHES "mp.*" AND NOT CXX_COMPILER_NAME MATCHES "CC" )
#  # message( FATAL_ERROR "MPI Compiler needed.\n"
#  #   "Please remove this build, export FC and CXX to a MPI compiler such as mpic++ or CC on Cray and rerun." )
#  set( HAVE_MPI OFF )
#endif()

### TRIANGULATION ...

ecbuild_add_option( FEATURE TRIANGULATION
                    DEFAULT OFF
                    DESCRIPTION "support for unstructured mesh generation"
                    REQUIRED_PACKAGES CGAL )

#   debug_var( CGAL_FOUND )
#   debug_var( CGAL_INCLUDE_DIRS )
#   debug_var( CGAL_LIBRARY )

### GRIB ...

#ecbuild_add_option( FEATURE GRIB
#					DEFAULT ON
#					DESCRIPTION "support for GRIB format encoding / decoding" )
#if( ENABLE_GRIB )
#	ecbuild_use_package( PROJECT grib_api VERSION 1.9.9 )
#	if( NOT GRIB_API_FOUND )
#		set( HAVE_GRIB 0 )
#	endif()
#endif()
#if( HAVE_GRIB )
#	add_definitions( -DHAVE_GRIB )
#	include_directories( ${GRIB_API_INCLUDE_DIRS} )
#	link_directories( ${GRIB_API_DIR}/lib )
#	list( APPEND TPL grib_api_f90 grib_api )
#else( )
#	message( STATUS "Grib not found : To enable, run \"cmake . -DGRIB_API_DIR=path/to/grib\"" )
#endif( )

### data dir ...

set( DATA_REPO  "ssh://git@software.ecmwf.int:7999/pan/data.git"
	 CACHE STRING "Location of the repository where large project files can be found")

message( STATUS "Retrieving meshes ..." )

include(GitExternal.cmake)
set( DATA ${CMAKE_BINARY_DIR}/data )
git_external( DIR ${DATA}
              REPO ${DATA_REPO}
			  TAG master )

add_definitions( -DATLAS_DATADIR="${DATA}")

### sources ...

add_subdirectory( src )

################################################################################################
# Symlink scripts

if( BUILD_SITE MATCHES "cc[abt].*" )
  message( STATUS "Creating symbolic link for submit_cray --> submit")
  execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
		${CMAKE_CURRENT_SOURCE_DIR}/tools/submit_cray ${CMAKE_CURRENT_BINARY_DIR}/submit )
endif( )

#######################################################################################
# Print Summary
#######################################################################################

ecbuild_install_project( NAME Atlas )

ecbuild_print_summary()

### final notes to user ...

message( STATUS "---------------------------------------------------------")
message( STATUS "NOTES:")
if( NOT CMAKE_BUILD_TYPE MATCHES "DEBUG" )
  message( STATUS " + Make sure the stack-size is big enough (ulimit -s unlimited)")
endif()
message( STATUS "---------------------------------------------------------")
