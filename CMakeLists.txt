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

#######################################################################################
# Configurations
#######################################################################################

ecbuild_find_python()

# OMP

ecbuild_add_option( FEATURE OMP
					DEFAULT ON
					DESCRIPTION "Enable/Disable OpenMP" )

debug_var( ENABLE_OMP )
debug_var( HAVE_OMP )

if( HAVE_OMP )
  add_definitions( -DHAVE_OMP )
endif( )

# MPI

ecbuild_add_option( FEATURE MPI
					DESCRIPTION "Enable/Disable MPI" )

# GRIB

ecbuild_add_option( FEATURE GRIB
					DESCRIPTION "support for GRIB format" )

if( ENABLE_GRIB )
	ecbuild_use_package( PROJECT grib_api VERSION 1.9.9 )
	if( NOT GRIB_API_FOUND )
		set( HAVE_GRIB 0 )
	endif()
endif()

if( HAVE_GRIB )
	add_definitions( -DHAVE_GRIB )
	include_directories( ${GRIB_API_INCLUDE_DIRS} )
	link_directories( ${GRIB_API_DIR}/lib )
	list( APPEND TPL grib_api_f90 grib_api )
else( )
	message( STATUS "Grib not found : To enable, run \"cmake . -DGRIB_API_DIR=path/to/grib\"" )
	set( GRIB OFF )
endif( )

### DATA

set( DATA_REPO     "ssh://git@software.ecmwf.int:7999/pan/data.git"
       CACHE STRING "Location of the repository where large project files can be found")

# Download external data from a git repository
message( STATUS "--------------------------------------------------------------------" )
message( STATUS "Retrieving meshes" )

include(GitExternal.cmake)
set( DATA ${CMAKE_BINARY_DIR}/data )
git_external( DIR ${DATA}
              REPO ${DATA_REPO}
              TAG master
             )

message( STATUS "--------------------------------------------------------------------" )
message( STATUS "Atlas configuration" )

get_filename_component(Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)
if( NOT Fortran_COMPILER_NAME MATCHES "mp.*" AND NOT Fortran_COMPILER_NAME MATCHES "ftn" )
  # message( FATAL_ERROR "MPI Compiler needed.\n"
  #   "Please remove this build, export FC and CXX to a MPI compiler such as mpif90 or ftn on Cray and rerun." )
  set( HAVE_MPI OFF )
endif( )
get_filename_component(CXX_COMPILER_NAME ${CMAKE_CXX_COMPILER} NAME)
if( NOT CXX_COMPILER_NAME MATCHES "mp.*" AND NOT CXX_COMPILER_NAME MATCHES "CC" )
  # message( FATAL_ERROR "MPI Compiler needed.\n"
  #   "Please remove this build, export FC and CXX to a MPI compiler such as mpic++ or CC on Cray and rerun." )
  set( HAVE_MPI OFF )
endif()

if( NOT CMAKE_BUILD_TYPE MATCHES "Debug" )
  add_definitions( -DNDEBUG )
endif( )

if( HAVE_MPI )
  add_definitions( -DHAVE_MPI )
endif( )

add_definitions( -DATLAS_DATADIR="${DATA}")

include_directories( src )

include( atlas_compiler_flags )

#######################################################################################
# Symlink scripts
#######################################################################################

if( BUILD_SITE MATCHES "cc[abt].*" )
  message( STATUS "Creating symbolic link for submit_cray --> submit")
  execute_process( COMMAND ln -sf ${CMAKE_CURRENT_SOURCE_DIR}/submit_cray ${CMAKE_CURRENT_BINARY_DIR}/submit )
endif( )

#######################################################################################
# Sources
#######################################################################################

add_subdirectory( src)

#######################################################################################
# Print Summary
#######################################################################################

ecbuild_install_project( NAME Atlas )

ecbuild_print_summary()

if ( NOT CMAKE_BUILD_TYPE MATCHES "DEBUG" )
  message( STATUS " + Make sure the stack-size is big enough (ulimit -s unlimited)")
endif()
message( STATUS "--------------------------------------------------------------------" )
