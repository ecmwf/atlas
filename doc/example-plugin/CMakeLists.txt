# (C) Copyright 2023- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
#
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.


cmake_minimum_required(VERSION 3.12 FATAL_ERROR)

find_package(ecbuild 3.4 REQUIRED HINTS ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/../ecbuild)

project(atlas-example-plugin LANGUAGES CXX)

find_package(atlas 0.25.0 REQUIRED)

atlas_create_plugin( atlas-example-plugin
    URL       https://atlas-example-plugin
    NAMESPACE int.ecmwf )

add_subdirectory(src)

ecbuild_install_project(NAME ${PROJECT_NAME})

