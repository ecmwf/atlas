# (C) Copyright 2013 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

ecbuild_add_executable(
    TARGET  atlas-parallel-interpolation
    SOURCES atlas-parallel-interpolation.cc
            PartitionedMesh.h
            PartitionedMesh.cc
    LIBS    atlas
    NOINSTALL
)

ecbuild_add_executable(
    TARGET  atlas-parallel-structured-interpolation
    SOURCES atlas-parallel-structured-interpolation.cc
    LIBS    atlas
    NOINSTALL
)

ecbuild_add_executable(
    TARGET  atlas-conservative-interpolation
    SOURCES atlas-conservative-interpolation.cc
    LIBS    atlas
    NOINSTALL
)

# Overcome Cray linking problem where the C++ library is not linked with the C library
list( APPEND NetCDF_CXX_EXTRA_LIBRARIES NetCDF::NetCDF_C )

ecbuild_add_option( FEATURE NETCDF
                    DESCRIPTION "Compile support the SCRIP format"
                    REQUIRED_PACKAGES "NetCDF COMPONENTS C CXX" )

ecbuild_add_executable(
    TARGET  atlas-global-matrix
    SOURCES atlas-global-matrix.cc
            ScripIO.h
            ScripIO.cc
    LIBS    atlas
)

target_compile_definitions( atlas-global-matrix PRIVATE ATLAS_HAVE_NETCDF=${HAVE_NETCDF} )
if( HAVE_NETCDF )
    target_link_libraries( atlas-global-matrix NetCDF::NetCDF_CXX )
endif()
