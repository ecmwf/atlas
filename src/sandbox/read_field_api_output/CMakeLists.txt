# (C) Copyright 2025- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

ecbuild_add_option( FEATURE HDF5
                    DESCRIPTION "HDF5 format"
                    REQUIRED_PACKAGES "HDF5 COMPONENTS C" )

if (HAVE_HDF5)
    ecbuild_add_executable(
        TARGET    atlas-read-field-api-output
        SOURCES   atlas-read-field-api-output.cc
        LIBS      atlas
                  ${HDF5_C_LIBRARIES}
        INCLUDES  ${HDF5_C_INCLUDE_DIRS}
        NOINSTALL
    )

    ecbuild_add_executable(
        TARGET    atlas-radintg
        SOURCES   atlas-radintg.cc
        LIBS      atlas
                  ${HDF5_C_LIBRARIES}
        INCLUDES  ${HDF5_C_INCLUDE_DIRS}
        NOINSTALL
    )

endif()
