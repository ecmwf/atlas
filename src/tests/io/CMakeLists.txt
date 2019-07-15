# (C) Copyright 2013 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

ecbuild_add_test( TARGET atlas_test_gmsh
  SOURCES  test_gmsh.cc ../TestMeshes.h
  LIBS     atlas
  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
)

ecbuild_add_test( TARGET atlas_test_pointcloud_io
  SOURCES   test_pointcloud_io.cc
  LIBS      atlas
  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
)

if( HAVE_FCTEST )

add_fctest( TARGET atlas_fctest_gmsh
  LINKER_LANGUAGE Fortran
  SOURCES         fctest_gmsh.F90
  LIBS            atlas_f
  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
)

endif()
