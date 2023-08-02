# (C) Copyright 2013 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

foreach(test
          test_bounding_box
          test_projection_LAEA
          test_rotation )

    ecbuild_add_test( TARGET atlas_${test} SOURCES ${test}.cc LIBS atlas ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT} )

endforeach()

ecbuild_add_test( TARGET atlas_test_jacobian SOURCES test_jacobian.cc LIBS atlas
    ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
                ATLAS_FPE=FE_OVERFLOW,FE_DIVBYZERO ) # FE_INVALID wrongly raised in release build of e.g. PGI 20.9

if(atlas_HAVE_PROJ)
    ecbuild_add_test( TARGET atlas_test_proj SOURCES test_proj.cc LIBS atlas ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT} )
endif()

if( HAVE_FCTEST )
    add_fctest( TARGET atlas_fctest_projection SOURCES fctest_projection.F90 LINKER_LANGUAGE Fortran LIBS atlas_f ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT} )
endif()


ecbuild_add_test( TARGET atlas_test_projection_variable_resolution SOURCES test_projection_variable_resolution.cc LIBS atlas ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
  CONDITION atlas_HAVE_ATLAS_FUNCTIONSPACE )

ecbuild_add_test( TARGET atlas_test_cubedsphere_projection
  MPI        6
  CONDITION  eckit_HAVE_MPI AND MPI_SLOTS GREATER_EQUAL 6 AND atlas_HAVE_ATLAS_FUNCTIONSPACE
  SOURCES test_cubedsphere_projection.cc
  LIBS atlas
  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
)
