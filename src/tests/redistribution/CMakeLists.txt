#
# (C) Crown Copyright 2021 Met Office
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
#

ecbuild_add_test( TARGET atlas_test_redistribution
  SOURCES   test_redistribution.cc
  MPI       8
  LIBS      atlas
  CONDITION eckit_HAVE_MPI AND MPI_SLOTS GREATER_EQUAL 8
  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
)
