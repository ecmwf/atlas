# (C) Copyright 1996-2014 ECMWF.

############################################################################################
# ATLAS
#
#		TODO:
#         * add renumbering to improve cache performance

cmake_minimum_required( VERSION 2.8.4 FATAL_ERROR )

project( atlas CXX )

set( CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/../ecbuild/cmake")
set( ECBUILD_DEFAULT_BUILD_TYPE Release )

include( ecbuild_system )

ecbuild_requires_macro_version( 1.3 )

ecbuild_declare_project()

################################################################################################
# options & dependencies

ecbuild_find_python()

### Fortran ...

ecbuild_add_option( FEATURE FORTRAN
					DEFAULT ON
					DESCRIPTION "provide Fortran bindings" )

if( ENABLE_FORTRAN )
	ecbuild_enable_fortran( REQUIRED MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/module )
endif()

### OMP ...

ecbuild_add_option( FEATURE OMP
					DEFAULT ON
					DESCRIPTION "support for OpenMP threaded parallelism" )

# !! TODO !! add test for OpenMP

### MPI ...

ecbuild_add_option( FEATURE MPI
					DEFAULT ON
					DESCRIPTION "support for MPI distributed parallelism" )

if( ENABLE_MPI )
   ecbuild_find_mpi()
   if( MPI_CXX_FOUND )

		add_definitions(${MPI_C_COMPILE_FLAGS})    # what is this ???
		include_directories(${MPI_C_INCLUDE_PATH})
		list( APPEND TPL ${MPI_C_LIBRARIES} )

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

### GRIB ...

ecbuild_add_option( FEATURE GRIB
					DEFAULT ON
					DESCRIPTION "support for GRIB format encoding / decoding" )
if( ENABLE_GRIB )
	ecbuild_use_package( PROJECT grib_api VERSION 1.9.9 )
	if( NOT GRIB_API_FOUND )
		set( HAVE_GRIB 0 )
	endif()
endif()
if( HAVE_GRIB )
	list( APPEND TPL grib_api_f90 grib_api )
endif( )

### data dir ...

set( DATA_REPO  "ssh://git@software.ecmwf.int:7999/pan/data.git"
	 CACHE STRING "Location of the repository where large project files can be found")

message( STATUS "Retrieving meshes ..." )

include(GitExternal.cmake)
set( ATLAS_DATADIR ${CMAKE_BINARY_DIR}/data )
git_external( DIR ${ATLAS_DATADIR}
              REPO ${DATA_REPO}
			  TAG master )

############################################################################################
# export package info


set( ATLAS_TPLS           ${TPL} )
set( ATLAS_INCLUDE_DIRS   ${CMAKE_CURRENT_SOURCE_DIR}/src ${CMAKE_CURRENT_BINARY_DIR}/src )
set( ATLAS_LIBRARIES      atlas )

if( ENABLE_FORTRAN )
	set( ATLAS_LIBRARIES  atlasf ${ATLAS_LIBRARIES} )
	list( APPEND ATLAS_INCLUDE_DIRS ${CMAKE_Fortran_MODULE_DIRECTORY} )
endif()

if( HAVE_CGAL )
	list( APPEND ATLAS_INCLUDE_DIRS ${CGAL_INCLUDE_DIRS} )
	list( APPEND ATLAS_LIBRARIES ${CGAL_LIBRARY} ${CGAL_3RD_PARTY_LIBRARIES} ${GMP_LIBRARIES} ${MPFR_LIBRARIES} )
endif()

get_directory_property( ATLAS_DEFINITIONS COMPILE_DEFINITIONS )

############################################################################################
# sources

include( atlas_compiler_flags )

include_directories( ${ATLAS_INCLUDE_DIRS} )

add_subdirectory( src )

############################################################################################
# symlink scripts

if( BUILD_SITE MATCHES "cc[abt].*" )
  message( STATUS "Creating symbolic link for submit_cray --> submit")
  execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
		${CMAKE_CURRENT_SOURCE_DIR}/tools/submit_cray ${CMAKE_BINARY_DIR}/bin/submit )
endif( )

############################################################################################
# finalize

set( ATLAS_INCLUDE_DIRS ${ATLAS_INCLUDE_DIRS} ${CMAKE_Fortran_MODULE_DIRECTORY} )

ecbuild_install_project( NAME Atlas )

ecbuild_print_summary()
