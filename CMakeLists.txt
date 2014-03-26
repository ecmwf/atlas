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
# TODO:
#	* add renumbering to improve cache performance
#	* add detection of OpenMP flags (mayeb move to ecbuild)

cmake_minimum_required( VERSION 2.8.4 FATAL_ERROR )

project( atlas CXX )

set( CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/../ecbuild/cmake")
set( ECBUILD_DEFAULT_BUILD_TYPE Release )

include( ecbuild_system )

ecbuild_requires_macro_version( 1.3 )

###############################################################################
# local project

set( PERSISTENT_NAMESPACE "eckit" CACHE INTERNAL "" )

ecbuild_use_package( PROJECT eckit VERSION 0.6 )

# init project

ecbuild_declare_project()

################################################################################################
# options & dependencies

ecbuild_find_python()

### Fortran ...

ecbuild_add_option( FEATURE FORTRAN_NUMBERING
					DEFAULT ON
					DESCRIPTION "data-structures are index with Fortran numbering" )

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
		cmake_add_c_flags( ${MPI_C_COMPILE_FLAGS} )
		cmake_add_cxx_flags( ${MPI_CXX_COMPILE_FLAGS} )
		include_directories(${MPI_C_INCLUDE_PATH})
	else()
		set( HAVE_MPI 0 )
   endif()
endif()

### grid

ecbuild_add_option( FEATURE GRID
					DEFAULT ON
					DESCRIPTION "support for C++ grid interface" )

if( ENABLE_GRID AND NOT ECKIT_FOUND )
	set( HAVE_GRID 0 )
	message( WARNING "atlas::grid requires eckit (not found)" )
else()
	ecbuild_add_cxx11_flags()
	ecbuild_check_cxx11( PRINT )
endif()

### tesselation

ecbuild_add_option( FEATURE TESSELATION
					DEFAULT ON
					DESCRIPTION "support for unstructured mesh generation"
					REQUIRED_PACKAGES CGAL )

### sandbox

ecbuild_add_option( FEATURE SANDBOX
					DEFAULT ON
					DESCRIPTION "build the sandbox stuff" )

############################################################################################
# export package info

set( ATLAS_INCLUDE_DIRS   ${CMAKE_CURRENT_SOURCE_DIR}/src ${CMAKE_CURRENT_BINARY_DIR}/src )
set( ATLAS_LIBRARIES      atlas )
set( ATLAS_TPLS           "" )

if( ENABLE_FORTRAN )
	list( APPEND ATLAS_INCLUDE_DIRS ${CMAKE_Fortran_MODULE_DIRECTORY} )
	set( ATLAS_LIBRARIES atlasf ${ATLAS_LIBRARIES} )
endif()

if( HAVE_TESSELATION AND CGAL_FOUND )
	list( APPEND CGAL_LIBRARIES ${CGAL_LIBRARY} ${CGAL_3RD_PARTY_LIBRARIES} ${GMP_LIBRARIES} ${MPFR_LIBRARIES} )
	list( APPEND ATLAS_TPLS CGAL )
endif()

if( HAVE_GRID )
	list( APPEND ATLAS_TPLS eckit )
	list( APPEND ATLAS_LIBRARIES atlas_grid )
endif()

get_directory_property( ATLAS_DEFINITIONS COMPILE_DEFINITIONS )

foreach( _tpl ${ATLAS_TPLS} )
	string( TOUPPER ${_tpl} TPL )
	list( APPEND ATLAS_EXTRA_DEFINITIONS   ${${TPL}_DEFINITIONS} )
	list( APPEND ATLAS_EXTRA_INCLUDE_DIRS  ${${TPL}_INCLUDE_DIRS} )
	list( APPEND ATLAS_EXTRA_LIBRARIES     ${${TPL}_LIBRARIES} )
endforeach()

############################################################################################
# sources

include( atlas_compiler_flags )

include_directories( ${ATLAS_INCLUDE_DIRS} ${ATLAS_EXTRA_INCLUDE_DIRS} )

add_subdirectory( src )

############################################################################################
# finalize

ecbuild_install_project( NAME Atlas )

ecbuild_print_summary()
