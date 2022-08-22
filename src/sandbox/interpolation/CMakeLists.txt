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
