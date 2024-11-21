# (C) Copyright 2024- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

############################################################################################
# hic

cmake_minimum_required( VERSION 3.18 FATAL_ERROR )

find_package( ecbuild 3.4 REQUIRED HINTS ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/../ecbuild ${CMAKE_CURRENT_SOURCE_DIR}/../../ecbuild )

################################################################################
# Initialise project

if( NOT DEFINED atlas_VERSION )
  set( atlas_VERSION 0.0.0 )
endif()

project( hic VERSION ${atlas_VERSION} LANGUAGES CXX )

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

################################################################################
# Features that can be enabled / disabled with -DENABLE_<FEATURE>

ecbuild_add_option( FEATURE CUDA DESCRIPTION "Use CUDA as backend to HIC" DEFAULT OFF )
ecbuild_add_option( FEATURE HIP  DESCRIPTION "Use HIP as backend to HIC"  DEFAULT OFF )

if( HAVE_CUDA )
  find_package(CUDAToolkit REQUIRED)
elseif( HAVE_HIP )
  find_package(hip CONFIG REQUIRED)
  find_package(hipsparse CONFIG REQUIRED)
endif()

add_subdirectory( src )
add_subdirectory( tests )

################################################################################
# Export and summarize

ecbuild_add_resources(
    TARGET hic-others
    SOURCES_PACK
        README.md
        CHANGELOG.md
        LICENSE
)

ecbuild_install_project( NAME hic )
ecbuild_print_summary()

