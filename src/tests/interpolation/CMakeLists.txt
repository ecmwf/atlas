# (C) Copyright 2013 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

ecbuild_add_test( TARGET atlas_test_Quad3D
  CONDITION ECKIT_HAVE_EIGEN
  SOURCES   test_Quad3D.cc
  LIBS      atlas
  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
)

ecbuild_add_test( TARGET atlas_test_interpolation_finite_element
  SOURCES   test_interpolation_finite_element.cc
  LIBS      atlas
  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
)

ecbuild_add_test( TARGET atlas_test_interpolation_cubic_prototype
  SOURCES  test_interpolation_cubic_prototype.cc CubicInterpolationPrototype.h
  LIBS     atlas
  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
)

ecbuild_add_executable( TARGET atlas_test_interpolation_structured2D
  SOURCES  test_interpolation_structured2D.cc
  LIBS     atlas
  NOINSTALL
)

if( ECBUILD-432-fixed )
# When ECBUILD-432 is fixed (already fixed in develop on top of ecbuild 3.0.0), we don't need TARGET_FILE
ecbuild_add_test( TARGET atlas_test_interpolation_bilinear
  COMMAND atlas_test_interpolation_structured2D ARGS --scheme linear
  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
)

ecbuild_add_test( TARGET atlas_test_interpolation_bicubic
  COMMAND atlas_test_interpolation_structured2D ARGS --scheme cubic
  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
)

ecbuild_add_test( TARGET atlas_test_interpolation_biquasicubic
  COMMAND atlas_test_interpolation_structured2D ARGS --scheme quasicubic
  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
)

else()
   # TO BE REMOVED!

  set( exe_atlas_test_interpolation_structured2D $<TARGET_FILE:atlas_test_interpolation_structured2D> )
  if( CMAKE_CROSSCOMPILING_EMULATOR )
    set( exe_atlas_test_interpolation_structured2D ${CMAKE_CROSSCOMPILING_EMULATOR} )
    set( arg_atlas_test_interpolation_structured2D $<TARGET_FILE:atlas_test_interpolation_structured2D> )
  endif()

  ecbuild_add_test( TARGET atlas_test_interpolation_bilinear
    COMMAND ${exe_atlas_test_interpolation_structured2D} ARGS ${arg_atlas_test_interpolation_structured2D} --scheme linear
    ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
  )

  ecbuild_add_test( TARGET atlas_test_interpolation_bicubic
    COMMAND ${exe_atlas_test_interpolation_structured2D} ARGS ${arg_atlas_test_interpolation_structured2D} --scheme cubic
    ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
  )

  ecbuild_add_test( TARGET atlas_test_interpolation_biquasicubic
    COMMAND ${exe_atlas_test_interpolation_structured2D} ARGS ${arg_atlas_test_interpolation_structured2D} --scheme quasicubic
    ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
  )
endif()



