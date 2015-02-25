
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
#  * add renumbering to improve cache performance
#  * add detection of OpenMP flags (mayeb move to ecbuild)

cmake_minimum_required( VERSION 2.8.4 FATAL_ERROR )

project( atlas CXX )

set( CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/../ecbuild/cmake")
set( ECBUILD_DEFAULT_BUILD_TYPE Release )

include( ecbuild_system )

ecbuild_requires_macro_version( 1.6 )

###############################################################################
# project

set( PERSISTENT_NAMESPACE "eckit" CACHE INTERNAL "" )

ecbuild_use_package( PROJECT eckit VERSION 0.8 REQUIRED )

foreach( _eckit_feature ${ECKIT_FEATURES} )
  set( ${_eckit_feature} 1 )
endforeach()

if( ECKIT_HAVE_GRIB )
   message( STATUS "eckit has GRIB support")
   set( ECKIT_GRIB_LIBRARIES eckit_grib )
endif()

if( ECKIT_HAVE_ECKIT_MPI )
   message( STATUS "eckit has mpi support")
endif()

if( ECKIT_HAVE_ECKIT_MPI )
   message( STATUS "eckit has mpistubs support")
endif()

ecbuild_declare_project()

# options & dependencies

ecbuild_find_python()

### Fortran ...

ecbuild_add_option( FEATURE FORTRAN
                    DESCRIPTION "provide Fortran bindings" )

if( ATLAS_HAVE_FORTRAN )

  ecbuild_enable_fortran( REQUIRED MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/module )

  ecbuild_add_option( FEATURE FCTEST
                      DEFAULT ${ENABLE_TESTS}
                      DESCRIPTION "Fortran Unit Test Framework"
                      REQUIRED_PACKAGES "PROJECT fctest VERSION 0.1" )

  set( Fortran Fortran )

endif()

### MPI ...

ecbuild_add_option( FEATURE MPI
                    DESCRIPTION "Support for MPI distributed parallelism"
                    REQUIRED_PACKAGES "MPI COMPONENTS CXX ${Fortran}" )

if( ATLAS_HAVE_MPI AND ECKIT_HAVE_ECKIT_MPI )
  message( STATUS "Using eckit_mpi" )
  set( ECKIT_MPI_LIBRARIES eckit_mpi )
elseif( NOT ATLAS_HAVE_MPI AND ECKIT_HAVE_ECKIT_MPISTUBS )
  message( STATUS "Using eckit_mpistubs" )
  set( ECKIT_MPI_LIBRARIES eckit_mpistubs )
  add_definitions( ECKIT_USE_MPISTUBS )
  if( ATLAS_HAVE_FORTRAN )
    ecbuild_use_package( PROJECT mpistubs )
    set( MPI_Fortran_LIBRARIES mpistubs_f )
    set( MPI_Fortran_INCLUDE_PATH ${MPISTUBS_INCLUDE_DIRS} )
  endif()
else()
  message( FATAL_ERROR "eckit_mpi or eckit_mpistubs not found" )
endif()


### OMP ...

ecbuild_add_option(
    FEATURE OMP
    DESCRIPTION "support for OpenMP threaded parallelism"
    REQUIRED_PACKAGES "OMP COMPONENTS CXX ${Fortran}" )

if( ATLAS_HAVE_FORTRAN )
  cmake_add_fortran_flags( "${OMP_Fortran_FLAGS}" )
endif()
cmake_add_cxx_flags( "${OMP_CXX_FLAGS}" )


### trans_api ...

ecbuild_add_option(
    FEATURE TRANS
    DESCRIPTION "Support for trasnformations between SH and Grid space"
    REQUIRED_PACKAGES "PROJECT trans_api VERSION 0.1 QUIET" )

### tesselation ...

set(Boost_USE_MULTITHREADED      ON )

ecbuild_add_option(
    FEATURE TESSELATION
    DESCRIPTION "support for unstructured mesh generation"
    REQUIRED_PACKAGES "CGAL QUIET" "Boost 1.45.0" )

if( HAVE_TESSELATION )
    list( APPEND CGAL_INCLUDE_DIRS ${Boost_INCLUDE_DIRS} )
    list( APPEND CGAL_LIBRARIES ${CGAL_LIBRARY} ${CGAL_3RD_PARTY_LIBRARIES} ${GMP_LIBRARIES} ${MPFR_LIBRARIES} ${Boost_THREAD_LIBRARY} ${Boost_SYSTEM_LIBRARY} )
endif()

### Type for Global indices and unique point id's

set( ATLAS_BITS_GLOBAL 64 )

### sandbox

ecbuild_add_option( FEATURE SANDBOX
                    DEFAULT OFF
                    DESCRIPTION "build the sandbox stuff" )

############################################################################################
# export package info

list( APPEND ATLAS_TPLS  eckit )

set( ATLAS_INCLUDE_DIRS   ${CMAKE_CURRENT_SOURCE_DIR}/src ${CMAKE_CURRENT_BINARY_DIR}/src )

list( APPEND ATLAS_LIBRARIES  atlas )

if( ATLAS_HAVE_FORTRAN )
  list( APPEND ATLAS_INCLUDE_DIRS ${CMAKE_Fortran_MODULE_DIRECTORY} )
  set( ATLAS_LIBRARIES atlas_f ${ATLAS_LIBRARIES} )
endif()

if( ATLAS_HAVE_TRANS )
  list( APPEND ATLAS_TPLS trans_api )
endif()

if( HAVE_TESSELATION )
  list( APPEND ATLAS_TPLS CGAL )
endif()

get_directory_property( ATLAS_DEFINITIONS COMPILE_DEFINITIONS )

foreach( _tpl ${ATLAS_TPLS} )
  string( TOUPPER ${_tpl} TPL )
  list( APPEND ATLAS_EXTRA_DEFINITIONS   ${${TPL}_DEFINITIONS}  ${${TPL}_TPL_DEFINITIONS}  )
  list( APPEND ATLAS_EXTRA_INCLUDE_DIRS  ${${TPL}_INCLUDE_DIRS} ${${TPL}_TPL_INCLUDE_DIRS} )
  list( APPEND ATLAS_EXTRA_LIBRARIES     ${${TPL}_LIBRARIES}    ${${TPL}_TPL_LIBRARIES}    )
endforeach()

############################################################################################
# sources

include( atlas_compiler_flags )

include_directories( ${ATLAS_INCLUDE_DIRS} ${ATLAS_EXTRA_INCLUDE_DIRS} )

add_subdirectory( src )

############################################################################################
# pkg-config

set( ATLAS_URL "https://software.ecmwf.int/wiki/display/ATLAS" )
set( ATLAS_DESCRIPTION "Atlas framework for parallel mesh datastructures" )

ecbuild_pkgconfig()

ecbuild_pkgconfig(
  NAME             atlas-c++
  LANGUAGES        CXX
  LIBRARIES        atlas
)

if( ATLAS_HAVE_FORTRAN )
  ecbuild_pkgconfig(
    NAME             atlas-fortran
    LANGUAGES        Fortran
    LIBRARIES        atlas_f
    NO_PRIVATE_INCLUDE_DIRS
  )
endif()

############################################################################################
# finalize

add_custom_target( atlas-readme SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/README.md )
ecbuild_install_project( NAME Atlas )

ecbuild_print_summary()
