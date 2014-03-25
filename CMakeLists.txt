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

ecbuild_use_package( PROJECT eckit    VERSION 0.6 REQUIRED )

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

### TRIANGULATION ...

ecbuild_add_option( FEATURE TRIANGULATION
                    DEFAULT OFF
                    DESCRIPTION "support for unstructured mesh generation"
                    REQUIRED_PACKAGES CGAL )

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
# finalize

set( ATLAS_INCLUDE_DIRS ${ATLAS_INCLUDE_DIRS} ${CMAKE_Fortran_MODULE_DIRECTORY} )

ecbuild_install_project( NAME Atlas )

ecbuild_print_summary()
