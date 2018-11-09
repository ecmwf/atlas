# (C) Copyright 2013 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.


if( HAVE_FCTEST )

  add_fctest( TARGET atlas_fctest_logging
    LINKER_LANGUAGE Fortran
    SOURCES         fctest_logging.F90
    LIBS            atlas_f
  )

  add_fctest( TARGET atlas_fctest_error
    LINKER_LANGUAGE Fortran
    SOURCES         fctest_error.F90
    LIBS            atlas_f
  )

 add_fctest( TARGET atlas_fctest_parametrisation
   CONDITION       NOT atlas_fctest_parametrisation_DISABLED
   LINKER_LANGUAGE Fortran
   SOURCES         fctest_parametrisation.F90
   LIBS            atlas_f
  )

  add_fctest( TARGET atlas_fctest_metadata
    LINKER_LANGUAGE Fortran
    SOURCES         fctest_metadata.F90
    LIBS            atlas_f
  )

endif()

foreach( test earth flags footprint indexview polygon point )
  ecbuild_add_test( TARGET atlas_test_${test}
    SOURCES test_${test}.cc
    LIBS atlas
  )
endforeach()

ecbuild_add_test( TARGET atlas_test_vector
  SOURCES  test_vector.cc
  LIBS     atlas
  CONDITION ATLAS_HAVE_GRIDTOOLS_STORAGE
)

ecbuild_add_test( TARGET atlas_test_metadata
  MPI        4
  CONDITION  ECKIT_HAVE_MPI
  SOURCES    test_metadata.cc
  LIBS       atlas
)

