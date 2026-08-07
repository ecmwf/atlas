# (C) Copyright 2013 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

ecbuild_add_test( TARGET atlas_test_algorithms_spectral
  MPI      4
  SOURCES   test_algorithms_spectral.cc
  CONDITION atlas_HAVE_ECTRANS AND eckit_HAVE_MPI AND ( transi_HAVE_MPI OR ectrans_HAVE_MPI ) AND MPI_SLOTS GREATER_EQUAL 4
  LIBS      atlas transi
  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
)
