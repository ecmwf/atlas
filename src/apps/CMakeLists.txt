# (C) Copyright 2013 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

ecbuild_add_executable(
  TARGET      atlas-main
  OUTPUT_NAME atlas
  SOURCES     atlas.cc
  LIBS        atlas )

ecbuild_add_executable(
  TARGET      atlas-meshgen
  SOURCES     atlas-meshgen.cc
  LIBS        atlas
  CONDITION   atlas_HAVE_ATLAS_FUNCTIONSPACE )

ecbuild_add_executable(
  TARGET      atlas-grids
  SOURCES     atlas-grids.cc
  LIBS        atlas
  CONDITION   atlas_HAVE_ATLAS_GRID )

ecbuild_add_executable(
  TARGET      atlas-grid-points
  SOURCES     atlas-grid-points.cc
  LIBS        atlas )

ecbuild_add_executable(
  TARGET      atlas-gaussian-latitudes
  SOURCES     atlas-gaussian-latitudes.cc
  LIBS        atlas
  CONDITION   atlas_HAVE_ATLAS_GRID )

ecbuild_add_executable(
  TARGET  atlas-interpolations
  SOURCES atlas-interpolations.cc
          ../atlas/output/ScripIO.h
          ../atlas/output/ScripIO.cc
  LIBS    atlas
)

# Overcome Cray linking problem where the C++ library is not linked with the C library
list( APPEND NetCDF_CXX_EXTRA_LIBRARIES NetCDF::NetCDF_C )

ecbuild_add_option( FEATURE NETCDF
                    DESCRIPTION "Compile support the SCRIP format"
                    REQUIRED_PACKAGES "NetCDF COMPONENTS C CXX" )
                    
target_compile_definitions( atlas-interpolations PRIVATE ATLAS_HAVE_NETCDF=${HAVE_NETCDF} )
if( HAVE_NETCDF )
  target_link_libraries( atlas-interpolations NetCDF::NetCDF_CXX )
endif()

add_custom_command( TARGET atlas-interpolations POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy
  ${CMAKE_CURRENT_SOURCE_DIR}/db_atlas_interpolations
  ${CMAKE_BINARY_DIR} )