# (C) Copyright 2013 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

############################################################################################
# ATLAS

cmake_minimum_required( VERSION 2.8.4 FATAL_ERROR )

set( CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/../ecbuild/cmake" )

set( ENABLE_OS_TYPES_TEST       OFF   CACHE BOOL " " )
set( ENABLE_OS_FUNCTIONS_TEST   OFF   CACHE BOOL " " )
set( ENABLE_LARGE_FILE_SUPPORT  OFF   CACHE BOOL " " )
set( ENABLE_OS_ENDINESS_TEST    OFF   CACHE BOOL " " )

project( atlas CXX )

include( ecbuild_system NO_POLICY_SCOPE )

ecbuild_requires_macro_version( 2.5 )

###############################################################################
# project

ecbuild_declare_project()

### eckit

ecbuild_use_package( PROJECT eckit VERSION 0.21.0 REQUIRED )
ecbuild_debug( "   ECKIT_FEATURES : [${ECKIT_FEATURES}]" )

# options & dependencies

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

ecbuild_find_python()

### Fortran ...

ecbuild_add_option( FEATURE FORTRAN
                    DESCRIPTION "Provide Fortran bindings"
                    REQUIRED_PACKAGES "PROJECT fckit VERSION 0.6.2" )

if( ATLAS_HAVE_FORTRAN )

  if( FCKIT_HAVE_ECKIT )

    ecbuild_enable_fortran( REQUIRED MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/module )
    set( Fortran Fortran )

    if( ENABLE_TESTS )
      set( HAVE_FCTEST ON )
      set( ATLAS_HAVE_FCTEST ON )
    else()
      set( HAVE_FCTEST OFF )
      set( ATLAS_HAVE_FCTEST OFF )
    endif()

  else()

    ecbuild_warn( "In order to compile atlas_f, fckit is required to be compiled with eckit. Turning off fortran." )
    set( HAVE_FORTRAN 0 )
    set( ATLAS_HAVE_FORTRAN 0 )
    unset( FCKIT_LIBRARIES )

  endif()
else()
  unset( FCKIT_LIBRARIES )
endif()

### MPI ...

# Atlas currently only required MPI for working around bug ECKIT-166 for MacOSX
if( ECKIT_HAVE_MPI AND EC_OS_NAME MATCHES "macosx" )
ecbuild_add_option( FEATURE MPI
                    DESCRIPTION "Support for MPI distributed parallelism"
                    REQUIRED_PACKAGES "MPI COMPONENTS CXX" )
endif()

if( NOT ECKIT_HAVE_MPI )
  ecbuild_warn("ecKit has been compiled without MPI. This causes Atlas to not be able to run parallel jobs.")
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

ecbuild_add_option( FEATURE FFTW
                    DESCRIPTION "Support for fftw"
                    REQUIRED_PACKAGES "FFTW COMPONENTS double" )

### trans ...

ecbuild_add_option( FEATURE TRANS
                    DESCRIPTION "Support for spectral transforms"
                    REQUIRED_PACKAGES "PROJECT transi VERSION 0.4.4 QUIET" )

### tesselation ...

set(Boost_USE_MULTITHREADED      ON )

ecbuild_add_option( FEATURE TESSELATION
                    DESCRIPTION "Support for unstructured mesh generation"
                    REQUIRED_PACKAGES "CGAL QUIET" "Boost 1.45.0" )

if( ATLAS_HAVE_TESSELATION )
    list( APPEND CGAL_INCLUDE_DIRS ${Boost_INCLUDE_DIRS} )
    if ( TARGET CGAL::CGAL )
      list( APPEND CGAL_LIBRARIES CGAL::CGAL ${CGAL_3RD_PARTY_LIBRARIES} ${GMP_LIBRARIES} ${MPFR_LIBRARIES} ${Boost_THREAD_LIBRARY} ${Boost_SYSTEM_LIBRARY} )
      # Reset INTERFACE_COMPILE_OPTIONS ( see ATLAS-193 )
      get_target_property( CGAL_COMPILE_FLAGS CGAL::CGAL INTERFACE_COMPILE_OPTIONS )
      set_target_properties( CGAL::CGAL PROPERTIES INTERFACE_COMPILE_OPTIONS "" )
    else()
      list( APPEND CGAL_LIBRARIES ${CGAL_LIBRARY} ${CGAL_3RD_PARTY_LIBRARIES} ${GMP_LIBRARIES} ${MPFR_LIBRARIES} ${Boost_THREAD_LIBRARY} ${Boost_SYSTEM_LIBRARY} )
    endif()
endif()

### GridTools storage module

ecbuild_add_option(
  FEATURE GRIDTOOLS_STORAGE
  DESCRIPTION "Arrays internally use GridTools storage layer"
  REQUIRED_PACKAGES "PROJECT gridtools_storage" )

if( ATLAS_HAVE_GRIDTOOLS_STORAGE )

  set( ATLAS_GRIDTOOLS_STORAGE_BACKEND_HOST 1 )
  set( ATLAS_GRIDTOOLS_STORAGE_BACKEND_CUDA 0 )

  if( GRIDTOOLS_HAVE_CUDA )
    ecbuild_info( "Gridtools found with CUDA support" )
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

if( ATLAS_HAVE_ACC )
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

ecbuild_add_option( FEATURE EIGEN
                    DESCRIPTION "Use Eigen linear algebra library"
                    REQUIRED_PACKAGES Eigen3 )

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
  if( NOT ATLAS_HAVE_INIT_SNAN )
    ecbuild_info( "Turning INIT_SNAN ON for Debug build" )
    set( ATLAS_HAVE_INIT_SNAN 1 )
  endif()
  if( NOT ATLAS_HAVE_BOUNDSCHECKING )
    ecbuild_info( "Turning BOUNDSCHECKING ON for Debug build" )
    set( ATLAS_HAVE_BOUNDSCHECKING 1 )
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

set( MIR_NEEDS_TRANSI TRUE CACHE BOOL "Deprecated" INTERNAL )

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

if( ATLAS_HAVE_PKGCONFIG )
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
        CHANGELOG.md
        LICENSE
)

if( ATLAS_HAVE_FORTRAN AND ECBUILD_INSTALL_FORTRAN_MODULES )
    install( DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/${CMAKE_CFG_INTDIR}
             DESTINATION module/atlas
             COMPONENT modules )
endif()

macro( filter_tpl _tpl )
  string( TOUPPER ${_tpl} _TPL )
  unset( ${_tpl}_INCLUDE_DIRS )
  unset( ${_TPL}_INCLUDE_DIRS )
  unset( ${_tpl}_INCLUDE_DIR )
  unset( ${_TPL}_INCLUDE_DIR )
  unset( ${_tpl}_LIBRARIES )
  unset( ${_TPL}_LIBRARIES )
  unset( ${_tpl}_LIBRARY )
  unset( ${_TPL}_LIBRARY )
  if( ${_tpl}_INCLUDE_DIRS )
    unset( ${_tpl}_INCLUDE_DIRS CACHE )
  endif()
  if( ${_TPL}_INCLUDE_DIRS )
    unset( ${_TPL}_INCLUDE_DIRS CACHE )
  endif()
  if( ${_tpl}_INCLUDE_DIR )
    unset( ${_tpl}_INCLUDE_DIR CACHE )
  endif()
  if( ${_TPL}_INCLUDE_DIR )
    unset( ${_TPL}_INCLUDE_DIR CACHE )
  endif()
  if( ${_tpl}_LIBRARIES )
    unset( ${_tpl}_LIBRARIES CACHE )
  endif()
  if( ${_TPL}_LIBRARIES )
    unset( ${_TPL}_LIBRARIES CACHE )
  endif()
  if( ${_tpl}_LIBRARY )
    unset( ${_tpl}_LIBRARY CACHE )
  endif()
  if( ${_TPL}_LIBRARY )
    unset( ${_TPL}_LIBRARY CACHE )
  endif()
endmacro()
filter_tpl( fckit )
filter_tpl( CGAL )
filter_tpl( FFTW )
filter_tpl( Eigen3 )
if( (NOT MIR_NEEDS_TRANSI) OR NOT( TRANSI_VERSION VERSION_LESS 0.6 ) )
  filter_tpl( transi )
endif()
get_target_property( eckit_intf_incl_dirs eckit INTERFACE_INCLUDE_DIRECTORIES )
set( eckit_can_be_filtered FALSE )
foreach( eckit_intf_incl_dir ${eckit_intf_incl_dirs} )
  string(FIND ${eckit_intf_incl_dir} "eckit" _found )
  if( NOT( _found EQUAL "-1" ) )
    set( eckit_can_be_filtered TRUE )
  endif()
endforeach()
if( eckit_can_be_filtered )
  filter_tpl(eckit)
endif()

ecbuild_install_project( NAME Atlas )

ecbuild_print_summary()

