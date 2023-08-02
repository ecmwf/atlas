# (C) Copyright 2023- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
#
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

configure_file(version.cc.in version.cc)

ecbuild_add_library(
    TARGET atlas-example-plugin
    TYPE SHARED
    SOURCES 
        Library.cc
        Library.h
        version.h
        ${CMAKE_CURRENT_BINARY_DIR}/version.cc
    HEADER_DESTINATION "include/atlas-example-plugin"
    PUBLIC_LIBS atlas
    PUBLIC_INCLUDES 
       $<INSTALL_INTERFACE:include/atlas-example-plugin>
       $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>
       $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/src> )

