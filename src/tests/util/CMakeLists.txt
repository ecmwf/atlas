# (C) Copyright 2013 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.


if( HAVE_FCTEST )

  add_fctest( TARGET atlas_fctest_functions
    LINKER_LANGUAGE Fortran
    SOURCES         fctest_functions.F90
    LIBS            atlas_f
    ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
  )

  add_fctest( TARGET atlas_fctest_logging
    LINKER_LANGUAGE Fortran
    SOURCES         fctest_logging.F90
    LIBS            atlas_f
    ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
  )

 add_fctest( TARGET atlas_fctest_parametrisation
   CONDITION       NOT atlas_fctest_parametrisation_DISABLED
   LINKER_LANGUAGE Fortran
   SOURCES         fctest_parametrisation.F90
   LIBS            atlas_f
   ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
  )

  add_fctest( TARGET atlas_fctest_metadata
    LINKER_LANGUAGE Fortran
    SOURCES         fctest_metadata.F90
    LIBS            atlas_f
    ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
  )

 add_fctest( TARGET atlas_fctest_geometry
    LINKER_LANGUAGE Fortran
    SOURCES         fctest_geometry.F90
    LIBS            atlas_f
    ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
  )

  add_fctest( TARGET atlas_fctest_kdtree
    LINKER_LANGUAGE Fortran
    SOURCES         fctest_kdtree.F90
    LIBS            atlas_f
    ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
  )
endif()

foreach( test util earth flags polygon point )
  ecbuild_add_test( TARGET atlas_test_${test}
    SOURCES test_${test}.cc
    LIBS atlas
    ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
  )
endforeach()

ecbuild_add_test( TARGET atlas_test_vector
  SOURCES  test_vector.cc
  LIBS     atlas
  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
)

ecbuild_add_test( TARGET atlas_test_metadata
  MPI        4
  CONDITION  eckit_HAVE_MPI AND MPI_SLOTS GREATER_EQUAL 4
  SOURCES    test_metadata.cc
  LIBS       atlas
  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
)

ecbuild_add_test( TARGET atlas_test_kdtree
  SOURCES  test_kdtree.cc
  LIBS     atlas
  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
  CONDITION atlas_HAVE_ATLAS_GRID
)

ecbuild_add_test( TARGET atlas_test_convexsphericalpolygon
  SOURCES test_convexsphericalpolygon.cc
  LIBS atlas
  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
)

ecbuild_add_test( TARGET atlas_test_unitsphere
  SOURCES test_unitsphere.cc
  LIBS atlas
  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
)

ecbuild_add_test( TARGET atlas_test_pack_vector_fields
  SOURCES test_pack_vector_fields.cc
  LIBS atlas
  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
)




