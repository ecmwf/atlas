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
# project

set( PERSISTENT_NAMESPACE "eckit" CACHE INTERNAL "" )

ecbuild_declare_project()

ecbuild_use_package( PROJECT eckit VERSION 0.6 )

# options & dependencies

ecbuild_find_python()

### Fortran ...

ecbuild_add_option( FEATURE FORTRAN_NUMBERING
                    DESCRIPTION "data-structures are index with Fortran numbering" )

ecbuild_add_option( FEATURE FORTRAN
                    DESCRIPTION "provide Fortran bindings" )

if( ENABLE_FORTRAN )
  ecbuild_enable_fortran( REQUIRED MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/module )

  ecbuild_add_option( FEATURE FCTEST
                      DEFAULT ${ENABLE_TESTS}
                      DESCRIPTION "Fortran Unit Test Framework"
                      REQUIRED_PACKAGES "PROJECT fctest VERSION 0.1.0" )
endif()

### OMP ...

ecbuild_add_option( FEATURE OMP
                    DESCRIPTION "support for OpenMP threaded parallelism" )

# !! TODO !! add test for OpenMP

### MPI ...

ecbuild_add_option( FEATURE MPI
                    DESCRIPTION "support for MPI distributed parallelism" )

if( ENABLE_MPI )
  ecbuild_find_mpi()
  if( MPI_CXX_FOUND )
    cmake_add_c_flags( " ${MPI_C_COMPILE_FLAGS}" )
    cmake_add_cxx_flags( " ${MPI_CXX_COMPILE_FLAGS}" )
    include_directories(${MPI_C_INCLUDE_PATH})
  else()
    set( HAVE_MPI 0 )
  endif()

  if( MPI_Fortran_FOUND )
    message( STATUS " MPI_Fortran_INCLUDE_PATH ${MPI_Fortran_INCLUDE_PATH}" )
    message( STATUS " MPI_Fortran_LIBRARIES ${MPI_Fortran_LIBRARIES}" )
    add_definitions(${MPI_Fortran_COMPILE_FLAGS})
    include_directories(${MPI_Fortran_INCLUDE_PATH})
  endif()
endif()

### grid

ecbuild_add_option( FEATURE GRID
                    DESCRIPTION "support for C++ grid interface" )

### tesselation

ecbuild_add_option( FEATURE TESSELATION
                    DESCRIPTION "support for unstructured mesh generation" )

if( ENABLE_TESSELATION )
  find_package( CGAL )

  if( CGAL_FOUND )
    list( APPEND CGAL_LIBRARIES ${CGAL_LIBRARY} ${CGAL_3RD_PARTY_LIBRARIES} ${GMP_LIBRARIES} ${MPFR_LIBRARIES} )
  endif()
   
  list( APPEND ATLAS_TPLS CGAL )
endif()

### Eigen

ecbuild_add_option( FEATURE EIGEN
                    DESCRIPTION "Eigen Linear Algebra Library"
                    REQUIRED_PACKAGES Eigen3 )

### Grib

ecbuild_add_option( FEATURE GRIB
                    DESCRIPTION "support for GRIB format encoding / decoding" )
debug_var( ENABLE_GRIB )
if( ENABLE_GRIB )
  ecbuild_use_package( PROJECT grib_api VERSION 1.9.9 )
  if( NOT GRIB_API_FOUND )
    set( HAVE_GRIB 0 )
  endif()
endif()
if( HAVE_GRIB )
  include_directories( ${GRIB_API_INCLUDE_DIRS} )
  list( APPEND TPL ${GRIB_API_LIBRARIES} )
endif( )

### sandbox

ecbuild_add_option( FEATURE SANDBOX
                    DESCRIPTION "build the sandbox stuff" )

############################################################################################
# export package info

list( APPEND ATLAS_TPLS eckit )
if( HAVE_GRIB )
  list( APPEND ATLAS_TPLS grib_api )
endif()

set( ATLAS_INCLUDE_DIRS   ${CMAKE_CURRENT_SOURCE_DIR}/src ${CMAKE_CURRENT_BINARY_DIR}/src )
set( ATLAS_LIBRARIES      atlas )

if( ENABLE_FORTRAN )
  list( APPEND ATLAS_INCLUDE_DIRS ${CMAKE_Fortran_MODULE_DIRECTORY} )
  set( ATLAS_LIBRARIES atlas_fortran atlas_mpl_fortran ${ATLAS_LIBRARIES} )
endif()

if( HAVE_GRID )
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
