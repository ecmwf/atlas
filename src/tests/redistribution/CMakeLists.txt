#
# (C) Crown Copyright 2021 Met Office
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
#

if( atlas_HAVE_ATLAS_FUNCTIONSPACE )

add_fctest( TARGET atlas_fctest_redistribution
  MPI 4
  CONDITION        eckit_HAVE_MPI
  LINKER_LANGUAGE  Fortran
  SOURCES          fctest_redistribution.F90
  LIBS             atlas_f
  ENVIRONMENT      ${ATLAS_TEST_ENVIRONMENT}
)

ecbuild_add_test( TARGET atlas_test_redistribution_structured
  SOURCES   test_redistribution_structured.cc
  MPI       8
  LIBS      atlas
  CONDITION eckit_HAVE_MPI AND MPI_SLOTS GREATER_EQUAL 8
  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
)

ecbuild_add_test( TARGET atlas_test_redistribution_generic
  SOURCES   test_redistribution_generic.cc
  MPI       8
  LIBS      atlas
  CONDITION eckit_HAVE_MPI AND MPI_SLOTS GREATER_EQUAL 8
  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
)

endif()
