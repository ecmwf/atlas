# (C) Copyright 1996-2014 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.


ecbuild_add_test( TARGET test_datastruct
  LINKER_LANGUAGE Fortran
  SOURCES  test_datastruct.F90
  LIBS     ftnunit datastruct
  COMMAND  ${CMAKE_SOURCE_DIR}/src/ftnunit/runtests.sh
  ARGS     ${CMAKE_CURRENT_BINARY_DIR}/test_datastruct
  CONDITION  OFF
)

ecbuild_add_test( TARGET test_gmsh
  SOURCES  test_gmsh.cpp
  LIBS     atlas
)

ecbuild_add_test( TARGET test_rgg
  SOURCES  test_rgg.cpp
  LIBS     atlas_meshgen atlas_testing
)

add_subdirectory( grid )
