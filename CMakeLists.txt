# (C) Copyright 2013 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

############################################################################################
# Atlas

cmake_minimum_required( VERSION 3.12 FATAL_ERROR )

find_package( ecbuild 3.4 REQUIRED HINTS ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/../ecbuild )

################################################################################
# Initialise project Atlas

project( atlas LANGUAGES CXX )

set( ATLAS_BITS_GLOBAL 64 CACHE STRING "bits used to define a global index (atlas::gidx_t)" )
set( ATLAS_BITS_LOCAL  32 CACHE STRING "bits used to define a local index (atlas::idx_t)" )

################################################################################
# Required packages

if( ATLAS_BITS_LOCAL MATCHES 64 )
  ecbuild_find_package( NAME eckit VERSION 1.17.1 REQUIRED )
else()
  ecbuild_find_package( NAME eckit VERSION 1.16 REQUIRED )
endif()

ecbuild_debug( "   eckit_FEATURES : [${eckit_FEATURES}]" )

add_subdirectory(atlas_io)
find_package(atlas_io)

################################################################################
# Features that can be enabled / disabled with -DENABLE_<FEATURE>

ecbuild_add_option( FEATURE ATLAS_GRID
                    DESCRIPTION "Build grid related features" )

ecbuild_add_option( FEATURE ATLAS_FIELD
                    DESCRIPTION "Build field and memory management related features" )

ecbuild_add_option( FEATURE ATLAS_FUNCTIONSPACE
                    DESCRIPTION "Build functionspace related features (mesh, functionspace, parallel)"
                    CONDITION atlas_HAVE_ATLAS_GRID AND atlas_HAVE_ATLAS_FIELD )

ecbuild_add_option( FEATURE ATLAS_INTERPOLATION
                    DESCRIPTION "Build interpolation related features"
                    CONDITION atlas_HAVE_ATLAS_FUNCTIONSPACE )

ecbuild_add_option( FEATURE ATLAS_TRANS
                    DESCRIPTION "Build transform related features"
                    CONDITION atlas_HAVE_ATLAS_FUNCTIONSPACE )

ecbuild_add_option( FEATURE ATLAS_NUMERICS
                    DESCRIPTION "Build numerics related features"
                    CONDITION atlas_HAVE_ATLAS_FUNCTIONSPACE )

ecbuild_add_option( FEATURE ECKIT_DEVELOP
                    DESCRIPTION "Used to enable new features or API depending on eckit develop branch, not yet in a tagged release"
                    DEFAULT OFF )


include( features/BOUNDSCHECKING )
include( features/FORTRAN )
include( features/CUDA )
include( features/MPI )
include( features/OMP )
include( features/FFTW )
include( features/ECTRANS )
include( features/TESSELATION )
include( features/GRIDTOOLS_STORAGE )
include( features/ACC )
include( features/EIGEN )
include( features/PROJ )
include( features/SANDBOX )
include( features/CLANG_TIDY )
include( features/INCLUDE_WHAT_YOU_USE )
include( features/INIT_SNAN )
include( features/DOCS )
include( features/ATLAS_RUN )

################################################################################
# sources

include( atlas_compile_flags )

add_subdirectory( hic )
add_subdirectory( src )
add_subdirectory( doc )

################################################################################
# Export and summarize

include( atlas_export )

ecbuild_print_summary()

