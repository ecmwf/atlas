# (C) Copyright 2013 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

############################################################################################
# Atlas

cmake_minimum_required( VERSION 3.6 FATAL_ERROR )

find_package( ecbuild 3.1.0 REQUIRED )

project( atlas LANGUAGES CXX )

### eckit

ecbuild_find_package( NAME eckit VERSION 1.4.0 REQUIRED )
ecbuild_debug( "   eckit_FEATURES : [${eckit_FEATURES}]" )

# options & dependencies

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

### Fortran ...

ecbuild_find_package( NAME fckit VERSION 0.6.2 COMPONENTS ECKIT )
ecbuild_add_option( FEATURE FORTRAN
                    DESCRIPTION "Provide Fortran bindings"
                    CONDITION fckit_FOUND )

if( atlas_HAVE_FORTRAN )

  if( fckit_HAVE_ECKIT )

    ecbuild_enable_fortran( REQUIRED MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/module )
    set( Fortran Fortran )

    if( HAVE_TESTS )
      set( HAVE_FCTEST ON )
      set( atlas_HAVE_FCTEST ON )
    else()
      set( HAVE_FCTEST OFF )
      set( atlas_HAVE_FCTEST OFF )
    endif()

  else()

    ecbuild_warn( "In order to compile atlas_f, fckit is required to be compiled with eckit. Turning off fortran." )
    set( HAVE_FORTRAN 0 )
    set( atlas_HAVE_FORTRAN 0 )

  endif()

endif()

ecbuild_find_python()

### MPI ...
if( NOT eckit_HAVE_MPI )
  ecbuild_warn("ecKit has been compiled without MPI. This causes Atlas to not be able to run parallel jobs.")
  set( atlas_HAVE_MPI 0 )
else()
  set( atlas_HAVE_MPI 1 )
endif()

### OMP ...
if( "${CMAKE_VERSION}" VERSION_LESS "3.11" )
  if( ENABLE_OMP )
    ecbuild_warn( "OpenMP only supported with CMake 3.11 onwards" )
  endif()
else()
  find_package( OpenMP COMPONENTS CXX ${Fortran} )
endif()
ecbuild_add_option( FEATURE OMP
                    DESCRIPTION "support for OpenMP shared memory parallelism"
                    CONDITION OpenMP_Fortran_FOUND OR OpenMP_CXX_FOUND )
ecbuild_add_option( FEATURE OMP_Fortran
                    DESCRIPTION "support for Fortran OpenMP shared memory parallelism"
                    CONDITION HAVE_OMP AND OpenMP_Fortran_FOUND )

ecbuild_add_option( FEATURE OMP_CXX
                    DESCRIPTION "support for CXX OpenMP shared memory parallelism"
                    CONDITION HAVE_OMP AND OpenMP_CXX_FOUND )
if( TARGET OpenMP::OpenMP_CXX )
  set( OMP_CXX OpenMP::OpenMP_CXX )
endif()
if( TARGET OpenMP::OpenMP_Fortran )
  set( OMP_Fortran OpenMP::OpenMP_Fortran )
endif()

### FFTW ...

find_package(FFTW COMPONENTS double QUIET )
ecbuild_add_option( FEATURE FFTW
                    DESCRIPTION "Support for fftw"
                    CONDITION FFTW_FOUND )

### trans ...

ecbuild_find_package( NAME transi VERSION 0.4.4 QUIET )
ecbuild_add_option( FEATURE TRANS
                    DESCRIPTION "Support for spectral transforms"
                    CONDITION transi_FOUND )

### tesselation ...

set(Boost_USE_MULTITHREADED      ON )

ecbuild_find_package( NAME CGAL QUIET)
ecbuild_find_package( NAME Boost VERSION 1.45.0 QUIET )
ecbuild_add_option( FEATURE TESSELATION
                    DESCRIPTION "Support for unstructured mesh generation"
                    CONDITION CGAL_FOUND AND Boost_FOUND )

if( atlas_HAVE_TESSELATION )
    list( APPEND CGAL_INCLUDE_DIRS ${Boost_INCLUDE_DIRS} )
    if ( TARGET CGAL::CGAL )
      list( APPEND CGAL_LIBRARIES CGAL::CGAL ${CGAL_3RD_PARTY_LIBRARIES} ${GMP_LIBRARIES} ${MPFR_LIBRARIES} ${Boost_THREAD_LIBRARY} ${Boost_SYSTEM_LIBRARY} )
      set( _cgal_target CGAL::CGAL )
      get_target_property(_aliased CGAL::CGAL ALIASED_TARGET)
      if(_aliased)
        set( _cgal_target ${_aliased} )
      endif()
      # Reset INTERFACE_COMPILE_OPTIONS ( see ATLAS-193 )
      get_target_property( CGAL_COMPILE_FLAGS ${_cgal_target} INTERFACE_COMPILE_OPTIONS )
      set_target_properties( ${_cgal_target} PROPERTIES INTERFACE_COMPILE_OPTIONS "" )
    else()
      list( APPEND CGAL_LIBRARIES ${CGAL_LIBRARY} ${CGAL_3RD_PARTY_LIBRARIES} ${GMP_LIBRARIES} ${MPFR_LIBRARIES} ${Boost_THREAD_LIBRARY} ${Boost_SYSTEM_LIBRARY} )
    endif()
endif()

### GridTools storage module

    ### GridTools may search for CUDA, which searches for "Threads"
    ### Set THREADS_HAVE_PTHREAD_ARG variable to false so that it can be recomputed based on
    ### THREADS_PREFER_PTHREAD_FLAG, in case other project had it on a different setting.
    ### This is certainly a CMake bug ( see ECKIT-426 )
    set( THREADS_HAVE_PTHREAD_ARG FALSE )
    if( NOT DEFINED THREADS_PREFER_PTHREAD_FLAG )
      set( THREADS_PREFER_PTHREAD_FLAG 1 )
    endif()

find_package( GridTools QUIET 
    HINTS ${GridTools_ROOT}/lib/cmake
          $ENV{GridTools_ROOT}/lib/cmake
          ${CMAKE_PREFIX_PATH}/lib/cmake
          ${CMAKE_INSTALL_PREFIX}/lib/cmake
          ${GridTools_BINARY_DIR} )
ecbuild_add_option(
  FEATURE GRIDTOOLS_STORAGE
  DESCRIPTION "Arrays internally use GridTools storage layer"
  CONDITION GridTools_FOUND )

if( atlas_HAVE_GRIDTOOLS_STORAGE )

  if( TARGET gridtools )
    # Fix gridtools export name. Otherwise "gridtools" will be in the INTERFACE_LINK_LIBRARIES of atlas in atlas-targets.cmake
    set_target_properties( gridtools PROPERTIES EXPORT_NAME GridTools::gridtools )
   endif()

  set( ATLAS_GRIDTOOLS_STORAGE_BACKEND_HOST 1 )
  set( ATLAS_GRIDTOOLS_STORAGE_BACKEND_CUDA 0 )

  if( GRIDTOOLS_HAS_BACKEND_CUDA )

    ecbuild_info( "GridTools found with CUDA support" )

    # Logic to check if we can use enable_language( CUDA )
    #   - CMake already supports it (as GridTools requires version that supports it)
    #   - ecbuild version 3.3 added support
    #   - overriding mechanism possible with cached "atlas_CUDA_LANGUAGE_ENABLED" variable
    if( DEFINED ecbuild_VERSION AND NOT ecbuild_VERSION VERSION_LESS 3.3 )
      set( atlas_CUDA_LANGUAGE_ENABLED_DEFAULT ON )
    else()
      set( atlas_CUDA_LANGUAGE_ENABLED_DEFAULT OFF )
    endif()
    set( atlas_CUDA_LANGUAGE_ENABLED ${atlas_CUDA_LANGUAGE_ENABLED_DEFAULT} CACHE BOOL "atlas enables CUDA language" )

    if( atlas_CUDA_LANGUAGE_ENABLED )
      enable_language( CUDA )
      ecbuild_info( "CUDA language enabled" )
    else()
      ecbuild_info("CUDA enabled through find_package(CUDA) instead of enable_language(CUDA)")
      find_package( CUDA )
    endif()

    set( ATLAS_GRIDTOOLS_STORAGE_BACKEND_HOST 0 )
    set( ATLAS_GRIDTOOLS_STORAGE_BACKEND_CUDA 1 )
  endif()

else()

  set( ATLAS_GRIDTOOLS_STORAGE_BACKEND_HOST 0 )
  set( ATLAS_GRIDTOOLS_STORAGE_BACKEND_CUDA 0 )

endif()


### OpenACC

set( ATLAS_ACC_CAPABLE FALSE )
if( ATLAS_GRIDTOOLS_STORAGE_BACKEND_CUDA )
  if( CMAKE_Fortran_COMPILER_ID MATCHES "PGI" )
    set( ATLAS_ACC_CAPABLE TRUE )
  endif()
endif()

ecbuild_add_option( FEATURE ACC
                    DESCRIPTION  "OpenACC capable data structures"
                    CONDITION ATLAS_ACC_CAPABLE )

if( atlas_HAVE_ACC )
  if( CMAKE_Fortran_COMPILER_ID MATCHES "PGI" )
    set( ACC_Fortran_FLAGS -acc -ta=tesla,nordc )
    set( ACC_C_FLAGS ${ACC_Fortran_FLAGS} )
    find_program( ACC_C_COMPILER NAMES pgcc HINTS ${PGI_DIR} ENV PGI_DIR PATH_SUFFIXES bin )
    if( NOT ACC_C_COMPILER )
      ecbuild_error( "Could not find OpenACC capable C compiler" )
    endif()
  endif()
endif()

### Eigen

ecbuild_find_package( NAME Eigen3 VERSION 3.3 QUIET )
ecbuild_add_option( FEATURE EIGEN
                    DESCRIPTION "Use Eigen linear algebra library"
                    CONDITION TARGET Eigen3::Eigen )

### Proj

ecbuild_find_package( NAME PROJ4 QUIET )
ecbuild_add_option( FEATURE PROJ
                    DESCRIPTION "PROJ-based projections"
                    DEFAULT OFF
                    CONDITION PROJ4_FOUND )

### Type for Global indices and unique point ids

set( ATLAS_BITS_GLOBAL 64 )
set( ATLAS_BITS_LOCAL 32 )

### Bounds checking
if( ${CMAKE_BUILD_TYPE} MATCHES "Debug" )
  set( DEFAULT_BOUNDSCHECKING ON )
  set( DEFAULT_INIT_SNAN ON )
else()
  set( DEFAULT_BOUNDSCHECKING OFF )
  set( DEFAULT_INIT_SNAN OFF )
endif()

ecbuild_add_option( FEATURE BOUNDSCHECKING
                    DEFAULT ${DEFAULT_BOUNDSCHECKING}
                    DESCRIPTION "Bounds checking for atlas::ArrayView and atlas::IndexView" )

ecbuild_add_option( FEATURE INIT_SNAN
                    DEFAULT ${DEFAULT_INIT_SNAN}
                    DESCRIPTION "Initialise atlas arrays with signaling_NaN (real types) or other invalid values (other types)" )

if( ${CMAKE_BUILD_TYPE} MATCHES "Debug" )
  if( NOT atlas_HAVE_INIT_SNAN )
    ecbuild_info( "Turning INIT_SNAN ON for Debug build" )
    set( atlas_HAVE_INIT_SNAN 1 )
  endif()
  if( NOT atlas_HAVE_BOUNDSCHECKING )
    ecbuild_info( "Turning BOUNDSCHECKING ON for Debug build" )
    set( atlas_HAVE_BOUNDSCHECKING 1 )
  endif()
endif()

### sandbox

ecbuild_add_option( FEATURE SANDBOX
                    DEFAULT OFF
                    DESCRIPTION "Build the sandbox stuff" )

### use of atlas-run for tests
ecbuild_add_option( FEATURE ATLAS_RUN
                    DEFAULT ON
                    DESCRIPTION "Use atlas/tools/atlas-run to run atlas tests" )

if( HAVE_ATLAS_RUN )
  set( MPIEXEC_EXECUTABLE ${CMAKE_CURRENT_SOURCE_DIR}/tools/atlas-run )
  set( MPIEXEC_NUMPROC_FLAG='-n' )
  set( CMAKE_CROSSCOMPILING_EMULATOR ${CMAKE_CURRENT_SOURCE_DIR}/tools/atlas-run )
endif()

################################################################################
# sources

include(CompileFlags)

include(FeatureClangTidy)

include(FeatureIncludeWhatYouUse)

add_subdirectory( src )

################################################################################
# export package info

if( TARGET atlas_f )
  list( APPEND ATLAS_LIBRARIES atlas_f )
endif()
list( APPEND ATLAS_LIBRARIES  atlas )

################################################################################
# pkg-config

ecbuild_add_option( FEATURE PKGCONFIG DESCRIPTION "Atlas pkgconfig" )
set( ATLAS_URL "https://software.ecmwf.int/wiki/display/ATLAS" )
set( ATLAS_DESCRIPTION "Atlas framework for parallel mesh datastructures" )

if( atlas_HAVE_PKGCONFIG )
  ecbuild_pkgconfig(
    NAME             atlas
    LIBRARIES        ${ATLAS_LIBRARIES}
  )

  ecbuild_pkgconfig(
    NAME             atlas-c++
    LANGUAGES        CXX
    LIBRARIES        atlas
  )

  if( atlas_HAVE_FORTRAN )
    ecbuild_pkgconfig(
      NAME             atlas-fortran
      LANGUAGES        Fortran
      LIBRARIES        atlas_f
      NO_PRIVATE_INCLUDE_DIRS
    )
  endif()
endif()

################################################################################
# documentation
if( ENABLE_DOCS ) 
  find_package(Latex)
endif()
ecbuild_add_option( FEATURE DOCS
                    DESCRIPTION "Atlas documentation"
                    DEFAULT OFF
                    CONDITION Latex_FOUND )

add_subdirectory( doc )

################################################################################
# finalize

ecbuild_add_resources(
    TARGET atlas-others
    SOURCES_PACK
        README.md
        CHANGELOG.md
        LICENSE
)

if( atlas_HAVE_FORTRAN AND ECBUILD_INSTALL_FORTRAN_MODULES )
    install( DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/${CMAKE_CFG_INTDIR}
             DESTINATION module/atlas
             COMPONENT modules )
endif()

set( atlas_REQUIRES_PRIVATE_DEPENDENCIES FALSE )
get_target_property( target_build_type atlas TYPE )
if( target_build_type STREQUAL STATIC_LIBRARY )
  set( atlas_REQUIRES_PRIVATE_DEPENDENCIES TRUE )
endif()

include( atlas_ecbuild2_compatibility )

ecbuild_install_project( NAME Atlas )

ecbuild_print_summary()

