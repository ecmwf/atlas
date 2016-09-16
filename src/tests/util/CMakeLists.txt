# (C) Copyright 1996-2016 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

set( FORTRAN_CALLBACK_LOGGING 0 )


if( HAVE_FCTEST )

  add_fctest( TARGET atlas_fctest_logging
    MPI             1
    LINKER_LANGUAGE Fortran
    SOURCES         fctest_logging.F90
    CONDITION ${FORTRAN_CALLBACK_LOGGING}
    LIBS            atlas_f
  )

  add_fctest( TARGET atlas_fctest_error
    MPI             1
    LINKER_LANGUAGE Fortran
    SOURCES         fctest_error.F90
    LIBS            atlas_f
  )

  add_fctest( TARGET atlas_fctest_parametrisation
    MPI             1
    LINKER_LANGUAGE Fortran
    SOURCES         fctest_parametrisation.F90
    LIBS            atlas_f
  )

  add_fctest( TARGET atlas_fctest_resource
    MPI             1
    LINKER_LANGUAGE Fortran
    SOURCES         fctest_resource.F90
    ARGS            -integer 10 -long 5000000000 -float 0.123456 -double 0.1234567890123456789 -string "hello world"
    LIBS            atlas_f
  )

endif()

ecbuild_add_test( TARGET atlas_test_flags
  BOOST
  MPI      1
  SOURCES  test_flags.cc
  LIBS     atlas
)

ecbuild_add_test( TARGET atlas_test_indexview
  BOOST
  MPI      1
  SOURCES  test_indexview.cc
  LIBS     atlas
)

ecbuild_add_test( TARGET atlas_test_metadata
  BOOST
  MPI      4
  SOURCES  test_metadata.cc
  LIBS     atlas
)
