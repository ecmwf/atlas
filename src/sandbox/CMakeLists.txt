# (C) Copyright 1996-2014 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

add_subdirectory( fortran_submodule )
add_subdirectory( fortran_modinc )
add_subdirectory( fortran_object )
add_subdirectory( example_fortran )

ecbuild_add_executable(
    TARGET  atlas-interpolation
    SOURCES atlas-interpolation.cc
    LIBS    atlas )

