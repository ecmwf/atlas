# (C) Copyright 1996-2016 ECMWF.
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
#  * add detection of OpenMP flags (maybe move to ecbuild)

cmake_minimum_required( VERSION 2.8.4 FATAL_ERROR )

set( CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/../ecbuild/cmake" )

set( ENABLE_OS_TYPES_TEST       OFF   CACHE BOOL " " )
set( ENABLE_OS_FUNCTIONS_TEST   OFF   CACHE BOOL " " )
set( ENABLE_LARGE_FILE_SUPPORT  OFF   CACHE BOOL " " )
set( ENABLE_OS_ENDINESS_TEST    OFF   CACHE BOOL " " )

project( atlas CXX )

file  ( STRINGS "VERSION.tex" ATLASVER )
set   ( ${PROJECT_NAME}_VERSION_STR  ${ATLASVER} )

include( ecbuild_system NO_POLICY_SCOPE )

ecbuild_requires_macro_version( 2.0 )

###############################################################################
# project

set( PERSISTENT_NAMESPACE "eckit" CACHE INTERNAL "" )

ecbuild_use_package( PROJECT eckit VERSION 0.12.3 REQUIRED )

ecbuild_debug( "   ECKIT_FEATURES : [${ECKIT_FEATURES}]" )

if( ECKIT_HAVE_ECKIT_MPI )
   ecbuild_debug( "eckit has mpi support" )
endif()

if( ECKIT_HAVE_ECKIT_MPISTUBS )
   ecbuild_debug( "eckit has mpistubs support" )
endif()

ecbuild_declare_project()

# options & dependencies

ecbuild_find_python()

### Fortran ...

ecbuild_add_option( FEATURE FORTRAN
                    DESCRIPTION "atlas: Provide Fortran bindings" )

if( ATLAS_HAVE_FORTRAN )

  ecbuild_enable_fortran( REQUIRED MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/module )

  ecbuild_add_option( FEATURE FCTEST
                      DEFAULT ${ENABLE_TESTS}
                      DESCRIPTION "atlas: Fortran Unit Test Framework"
                      REQUIRED_PACKAGES "PROJECT fckit VERSION 0.1" )

  set( Fortran Fortran )

endif()

### MPI ...

ecbuild_add_option( FEATURE MPI
                    DESCRIPTION "atlas: Support for MPI distributed parallelism"
                    REQUIRED_PACKAGES "MPI COMPONENTS CXX ${Fortran}" )

if( ATLAS_HAVE_MPI AND ECKIT_HAVE_ECKIT_MPI )
  ecbuild_debug( "Using eckit_mpi" )
  set( ECKIT_MPI_LIBRARIES eckit_mpi )
elseif( NOT ATLAS_HAVE_MPI AND ECKIT_HAVE_ECKIT_MPISTUBS )
  ecbuild_debug( "Using eckit_mpistubs" )
  set( ECKIT_MPI_LIBRARIES eckit_mpistubs )
  add_definitions( "-DECKIT_USE_MPISTUBS" )
  if( ATLAS_HAVE_FORTRAN )
    ecbuild_use_package( PROJECT mpistubs REQUIRED )
    set( MPI_Fortran_LIBRARIES mpistubs_f )
    set( MPI_Fortran_INCLUDE_PATH ${MPISTUBS_INCLUDE_DIRS} )
  endif()
elseif( NOT ATLAS_HAVE_MPI AND NOT ECKIT_HAVE_ECKIT_MPISTUBS )
  ecbuild_error("ATLAS_HAVE_MPI=OFF, but cannot build since eckit_mpistubs not found")
else()
  ecbuild_info( "ATLAS_HAVE_MPI            : ${ATLAS_HAVE_MPI}" )
  ecbuild_info( "ECKIT_HAVE_ECKIT_MPI      : ${ECKIT_HAVE_ECKIT_MPI} " )
  ecbuild_info( "ECKIT_HAVE_ECKIT_MPISTUBS : ${ECKIT_HAVE_ECKIT_MPISTUBS} " )
  ecbuild_error( "Invalid combination of above variables" )
endif()
if( ATLAS_HAVE_FORTRAN )
  ecbuild_debug( "MPI_Fortran_INCLUDE_PATH : ${MPI_Fortran_INCLUDE_PATH}" )
  ecbuild_debug( "MPI_Fortran_LIBRARIES    : ${MPI_Fortran_LIBRARIES}" )
endif()


### OMP ...

ecbuild_add_option(
    FEATURE OMP
    DESCRIPTION "atlas: Support for OpenMP threaded parallelism"
    REQUIRED_PACKAGES "OMP COMPONENTS CXX ${Fortran}" )

if( HAVE_OMP )
  ecbuild_enable_omp()
else()
  ecbuild_enable_ompstubs()
endif()

### trans ...

ecbuild_add_option(
    FEATURE TRANS
    DESCRIPTION "atlas: Support for spectral transforms"
    REQUIRED_PACKAGES "PROJECT transi VERSION 0.3.2 QUIET" )

### tesselation ...

set(Boost_USE_MULTITHREADED      ON )

ecbuild_add_option(
    FEATURE TESSELATION
    DESCRIPTION "atlas: Support for unstructured mesh generation"
    REQUIRED_PACKAGES "CGAL QUIET" "Boost 1.45.0" )

if( HAVE_TESSELATION )
    list( APPEND CGAL_INCLUDE_DIRS ${Boost_INCLUDE_DIRS} )
    if ( TARGET CGAL::CGAL )
      list( APPEND CGAL_LIBRARIES CGAL::CGAL ${CGAL_3RD_PARTY_LIBRARIES} ${GMP_LIBRARIES} ${MPFR_LIBRARIES} ${Boost_THREAD_LIBRARY} ${Boost_SYSTEM_LIBRARY} )
    else()
      list( APPEND CGAL_LIBRARIES ${CGAL_LIBRARY} ${CGAL_3RD_PARTY_LIBRARIES} ${GMP_LIBRARIES} ${MPFR_LIBRARIES} ${Boost_THREAD_LIBRARY} ${Boost_SYSTEM_LIBRARY} )
    endif()
endif()

### FIND GRIDTOOLS STORAGE
ecbuild_add_option(
    FEATURE GRIDTOOLS_STORAGE
    DESCRIPTION "atlas: Arrays internally use GridTools storage layer"
    REQUIRED_PACKAGES "gridtools_storage" )

ecbuild_add_cxx11_flags()

if( ATLAS_HAVE_GRIDTOOLS_STORAGE )
    set( ENABLE_GPU "OFF" CACHE BOOL "Compile with GPU support (CUDA)" )

    if(ENABLE_GPU)
        find_package(CUDA REQUIRED)
        string(REPLACE "." "" CUDA_VERSION ${CUDA_VERSION})
        add_definitions(-DCUDA_VERSION=${CUDA_VERSION})
        add_definitions(-DENABLE_GPU)
        set(CUDA_PROPAGATE_HOST_FLAGS OFF)
        if( ${CUDA_VERSION} VERSION_GREATER "60")
            set(CUDA_NVCC_FLAGS "--std=c++11" "--relaxed-constexpr" "${CUDA_NVCC_FLAGS}")
        else()
            message(FATAL_ERROR "CUDA 6.0 or lower does not support C++11 (disabling)")
        endif()

        set( CUDA_ARCH "sm_35" CACHE STRING "Compute capability for CUDA" )
        set( CUDA_SEPARABLE_COMPILATION ON)

        include_directories(SYSTEM ${CUDA_INCLUDE_DIRS})

        set(exe_LIBS "${CUDA_CUDART_LIBRARY}" "${exe_LIBS}" )
        set(CUDA_SEPARABLE_COMPILATION ON)
        # adding the additional nvcc flags
        set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS}" "-arch=${CUDA_ARCH}" "-Xcudafe" "--diag_suppress=dupl_calling_convention")
        set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS}" "-Xcudafe" "--diag_suppress=code_is_unreachable" "-Xcudafe")
        set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS}" "--diag_suppress=implicit_return_from_non_void_function" "-Xcudafe")
        set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS}" "--diag_suppress=calling_convention_not_allowed" "-Xcudafe")
        set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS}" "--diag_suppress=conflicting_calling_conventions")

    endif()
endif()

### Eigen

ecbuild_add_option( FEATURE EIGEN
                    DESCRIPTION "Eigen support"
                    TYPE RECOMMENDED
                    PURPOSE "Eigen backend: accelerated small matrix-vector operations"
                    REQUIRED_PACKAGES Eigen3 )

### Type for Global indices and unique point id's

set( ATLAS_BITS_GLOBAL 64 )

### Bounds checking
if( ${CMAKE_BUILD_TYPE} MATCHES "Debug" )
  set( DEFAULT_BOUNDSCHECKING ON )
else()
  set( DEFAULT_BOUNDSCHECKING OFF )
endif()
ecbuild_add_option( FEATURE BOUNDSCHECKING
                    DEFAULT ${DEFAULT_BOUNDSCHECKING}
                    DESCRIPTION "atlas: Bounds checking for atlas::ArrayView and atlas::IndexView" )
if( ${CMAKE_BUILD_TYPE} MATCHES "Debug" AND NOT ATLAS_HAVE_BOUNDSCHECKING )
  ecbuild_info( "Turning BOUNDSCHECKING ON for Debug build" )
  set( ATLAS_HAVE_BOUNDSCHECKING 1 )
endif()

### sandbox

ecbuild_add_option( FEATURE SANDBOX
                    DEFAULT OFF
                    DESCRIPTION "Build the sandbox stuff" )

################################################################################
# export package info

set( ATLAS_INCLUDE_DIRS   ${CMAKE_CURRENT_SOURCE_DIR}/src ${CMAKE_CURRENT_BINARY_DIR}/src )

list( APPEND ATLAS_LIBRARIES  atlas )

if( ATLAS_HAVE_FORTRAN )
  list( APPEND ATLAS_INCLUDE_DIRS ${CMAKE_Fortran_MODULE_DIRECTORY} )
  set( ATLAS_LIBRARIES atlas_f ${ATLAS_LIBRARIES} )
endif()

get_directory_property( ATLAS_DEFINITIONS COMPILE_DEFINITIONS )

foreach( _tpl ${ATLAS_TPLS} )
  string( TOUPPER ${_tpl} TPL )
  list( APPEND ATLAS_EXTRA_DEFINITIONS   ${${TPL}_DEFINITIONS}  ${${TPL}_TPL_DEFINITIONS}  )
  list( APPEND ATLAS_EXTRA_INCLUDE_DIRS  ${${TPL}_INCLUDE_DIRS} ${${TPL}_TPL_INCLUDE_DIRS} )
  list( APPEND ATLAS_EXTRA_LIBRARIES     ${${TPL}_LIBRARIES}    ${${TPL}_TPL_LIBRARIES}    )
endforeach()

ecbuild_debug_var( ATLAS_EXTRA_INCLUDE_DIRS )
ecbuild_debug_var( ATLAS_EXTRA_LIBRARIES )


################################################################################
# sources

include_directories( ${ATLAS_INCLUDE_DIRS} ${ATLAS_EXTRA_INCLUDE_DIRS} )
add_subdirectory( src )

################################################################################
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

################################################################################
# documentation
ecbuild_add_option( FEATURE DOCS
                    DESCRIPTION "Atlas documentation"
                    DEFAULT OFF
                    REQUIRED_PACKAGES "Latex" )

if(ATLAS_HAVE_DOCS)
    add_subdirectory( doc )
endif()

################################################################################
# finalize

ecbuild_add_resources(
    TARGET atlas-others
    SOURCES_PACK
        README.md
        LICENSE
)

ecbuild_install_project( NAME Atlas )

if( ATLAS_HAVE_FORTRAN )
    install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/${CMAKE_CFG_INTDIR} DESTINATION ${INSTALL_INCLUDE_DIR} )
endif()

ecbuild_print_summary()

