# (C) Copyright 2013 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

ecbuild_add_test( TARGET atlas_test_haloexchange_adjoint
  MPI        3
  CONDITION  eckit_HAVE_MPI
  SOURCES    test_haloexchange_adjoint.cc
  LIBS       atlas
  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
)

ecbuild_add_test( TARGET atlas_test_setcomm
  MPI        3
  CONDITION  eckit_HAVE_MPI
  SOURCES    test_setcomm.cc
  LIBS       atlas
  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
)

ecbuild_add_test( TARGET atlas_test_haloexchange
  MPI        3
  CONDITION  eckit_HAVE_MPI
  SOURCES    test_haloexchange.cc
  LIBS       atlas
  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
)

ecbuild_add_test( TARGET atlas_test_haloexchange_on_device
  MPI        3
  CONDITION  eckit_HAVE_MPI AND atlas_HAVE_GPU
  SOURCES    test_haloexchange.cc
  LIBS       atlas
  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT} ATLAS_RUN_NGPUS=1
)
if( TEST atlas_test_haloexchange_on_device )
  set_tests_properties ( atlas_test_haloexchange_on_device PROPERTIES LABELS "gpu")
endif()

ecbuild_add_test( TARGET atlas_test_gather
  MPI        3
  CONDITION  eckit_HAVE_MPI
  SOURCES    test_gather.cc
  LIBS       atlas
  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
)

ecbuild_add_test( TARGET atlas_test_omp_sort
  OMP        8
  SOURCES    test_omp_sort.cc
  LIBS       atlas
  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
)

ecbuild_add_test( TARGET atlas_test_omp_copy
  OMP        8
  SOURCES    test_omp_copy.cc
  LIBS       atlas
  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
)

