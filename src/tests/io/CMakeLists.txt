# (C) Copyright 1996-2016 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

ecbuild_add_test( TARGET atlas_test_gmsh
  BOOST
  MPI      1
  SOURCES  test_gmsh.cc ../TestMeshes.h
  LIBS     atlas
)

ecbuild_add_test( TARGET atlas_test_pointcloud
  BOOST
  MPI       1
  SOURCES   test_pointcloud.cc
  LIBS      atlas
)
