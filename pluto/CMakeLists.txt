# (C) Copyright 2024- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

cmake_minimum_required( VERSION 3.18 FATAL_ERROR )

find_package( ecbuild 3.8 REQUIRED HINTS ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/../ecbuild ${CMAKE_CURRENT_SOURCE_DIR}/../../ecbuild )

if( NOT DEFINED atlas_VERSION )
  set( atlas_VERSION 0.0.0 )
endif()

project( pluto VERSION ${atlas_VERSION} LANGUAGES CXX )

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(PLUTO_HAVE_PMR 0)
### Normally we should be auto-detecting std::pmr (c++17) for polymorphic memory resource
#   However, CUDA or HIP compiler seems to not find the right headers,
#   so for now we don't rely on it
#
#    include(CheckIncludeFileCXX)
#    check_include_file_cxx(memory_resource PLUTO_HAVE_PMR)
#    if( PLUTO_HAVE_PMR )
#      set( PLUTO_HAVE_PMR 1 )
#    else()
#      set( PLUTO_HAVE_PMR 0 )
#    endif()

find_package(hic QUIET)
if( NOT hic_FOUND )
  add_subdirectory(../hic ${CMAKE_CURRENT_BINARY_DIR}/hic )
  find_package(hic REQUIRED)
endif()

set( HAVE_CUDA ${hic_HAVE_CUDA} )
set( HAVE_HIP  ${hic_HAVE_HIP}  )

include("cmake/pluto_hic_macros.cmake")

if( HAVE_CUDA )
  ecbuild_info("pluto has GPU support via CUDA")
  set( CMAKE_CUDA_RESOLVE_DEVICE_SYMBOLS ON )
  set( CMAKE_CUDA_SEPARABLE_COMPILATION  ON )
  enable_language(CUDA)
elseif( HAVE_HIP )
  ecbuild_info("pluto has GPU support via HIP")
  set(CMAKE_HIP_STANDARD 17)
  set(CMAKE_HIP_STANDARD_REQUIRED ON)
  enable_language(HIP)
else()
  ecbuild_info("pluto has no GPU support")
endif()

ecbuild_add_option( FEATURE WARNINGS  DESCRIPTION "Add warnings to compiler"   DEFAULT ON  )

if( HAVE_CUDA )
  find_package(CUDAToolkit REQUIRED)
elseif( HAVE_HIP )
  find_package(hip       CONFIG REQUIRED)
  find_package(hipsparse CONFIG REQUIRED)
endif()

# activate warnings, ecbuild macros check the compiler recognises the options
if(HAVE_WARNINGS)
  ecbuild_add_cxx_flags("-Wall"      NO_FAIL NAME pluto_cxx_Wall QUIET)
  ecbuild_add_cxx_flags("-Wextra"    NO_FAIL NAME pluto_cxx_Wextra QUIET)
  ecbuild_add_cxx_flags("-Wpedantic" NO_FAIL NAME pluto_cxx_Wpedantic QUIET)
endif()

ecbuild_add_option( FEATURE WARNING_AS_ERROR
                    DEFAULT OFF
                    DESCRIPTION "Treat compile warning as error" )

if(HAVE_WARNING_AS_ERROR)
  ecbuild_add_cxx_flags("-Werror" NO_FAIL NAME pluto_cxx_warning_as_error QUIET)
endif()


ecbuild_add_option( FEATURE FORTRAN DESCRIPTION "Fortran API for pluto" )
if( HAVE_FORTRAN )
    ecbuild_enable_fortran( REQUIRED MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/module )
endif()

add_subdirectory(src)
add_subdirectory(examples)
add_subdirectory(tests)

ecbuild_install_project( NAME pluto )
ecbuild_print_summary()
