# (C) Copyright 2013 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

if( HAVE_FCTEST )

  if( ATLAS_HAVE_TRANS )
    add_fctest( TARGET atlas_fctest_trans
      LINKER_LANGUAGE Fortran
      SOURCES         fctest_trans.F90
      LIBS            atlas_f
      CONDITION       ECKIT_HAVE_MPI AND TRANSI_HAVE_MPI
      MPI             4
      ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
    )
  endif()

  if( ATLAS_HAVE_TRANS )
    add_fctest( TARGET atlas_fctest_trans_invtrans_grad
      LINKER_LANGUAGE Fortran
      SOURCES         fctest_trans_invtrans_grad.F90
      LIBS            atlas_f
      ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
    )
  endif()

  add_fctest( TARGET atlas_fctest_trans_unstructured
    LINKER_LANGUAGE Fortran
    SOURCES         fctest_trans_unstructured.F90
    LIBS            atlas_f
  )

endif()

ecbuild_add_test( TARGET atlas_test_trans
  MPI      4
  SOURCES   test_trans.cc
  CONDITION ATLAS_HAVE_TRANS AND ECKIT_HAVE_MPI AND TRANSI_HAVE_MPI
  LIBS      atlas transi
  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
)

ecbuild_add_test( TARGET atlas_test_trans_serial
  SOURCES   test_trans.cc
  CONDITION ATLAS_HAVE_TRANS
  LIBS      atlas transi
  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
)

ecbuild_add_test( TARGET atlas_test_trans_invtrans_grad
  SOURCES   test_trans_invtrans_grad.cc
  CONDITION ATLAS_HAVE_TRANS
  LIBS      atlas transi
  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
)

# Note: This duplication is needed because transi is a private library of atlas,
#       and this tests needs access to the transi include directories.
# ToDo: Fix this inside the test code so that we don't directly need to include transi headers.
if( ATLAS_HAVE_TRANS )
ecbuild_add_test( TARGET atlas_test_transgeneral
  SOURCES   test_transgeneral.cc
  LIBS      atlas transi
  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT} ATLAS_TRACE_REPORT=1
  CONDITION ATLAS_HAVE_FFTW
)
else()
ecbuild_add_test( TARGET atlas_test_transgeneral
  SOURCES   test_transgeneral.cc
  LIBS      atlas
  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT} ATLAS_TRACE_REPORT=1
  CONDITION ATLAS_HAVE_FFTW
)
endif()


ecbuild_add_test( TARGET atlas_test_trans_localcache
  SOURCES   test_trans_localcache.cc
  LIBS      atlas
  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT} ATLAS_TRACE_REPORT=1
)

