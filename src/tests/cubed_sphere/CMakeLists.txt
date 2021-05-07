#
# (C) Crown Copyright 2021 Met Office
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
#

ecbuild_add_test( TARGET atlas_test_cubed_sphere
  SOURCES   test_cubed_sphere.cc
  LIBS      atlas
  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
)
