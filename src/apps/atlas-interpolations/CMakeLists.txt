# (C) Copyright 2025 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.


ecbuild_add_executable(
  TARGET  atlas-interpolations
  SOURCES atlas-interpolations.cc
          ScripIO.h
          ScripIO.cc
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

