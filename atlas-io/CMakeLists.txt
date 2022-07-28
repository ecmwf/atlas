# (C) Copyright 2021 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

############################################################################################
# Atlas-IO

cmake_minimum_required( VERSION 3.12 FATAL_ERROR )

find_package( ecbuild 3.4 REQUIRED HINTS ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/../ecbuild )

################################################################################
# Initialise project Atlas

project( atlas-io VERSION ${atlas_VERSION} LANGUAGES CXX )

################################################################################
# Required packages

ecbuild_find_package( NAME eckit REQUIRED )

ecbuild_debug( "   eckit_FEATURES : [${eckit_FEATURES}]" )

################################################################################
# Features that can be enabled / disabled with -DENABLE_<FEATURE>

ecbuild_add_option( FEATURE CXX17
                    DESCRIPTION "Use C++17 standard"
                    DEFAULT OFF )

################################################################################
# sources

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

check_cxx_source_compiles( "#include <cxxabi.h>\n int main() { char * type; int status; char * r = abi::__cxa_demangle(type, 0, 0, &status); }"
    ATLAS_IO_HAVE_CXXABI_H )

test_big_endian( _BIG_ENDIAN )

if( _BIG_ENDIAN )
  set( ATLAS_IO_BIG_ENDIAN    1 )
  set( ATLAS_IO_LITTLE_ENDIAN 0 )
else()
  set( ATLAS_IO_BIG_ENDIAN    0 )
  set( ATLAS_IO_LITTLE_ENDIAN 1 )
endif()


add_subdirectory( src )
add_subdirectory( tests )

################################################################################
# Export and summarize

ecbuild_add_resources(
    TARGET atlas-io-others
    SOURCES_PACK
        README.md
        CHANGELOG.md
        LICENSE
)

ecbuild_install_project( NAME Atlas-IO )
ecbuild_print_summary()

