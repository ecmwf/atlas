# (C) Copyright 2013 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

ecbuild_add_test( TARGET atlas_test_Quad3D
  CONDITION eckit_EIGEN_FOUND
  SOURCES   test_Quad3D.cc
  LIBS      atlas
  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
)

ecbuild_add_test( TARGET atlas_test_Quad2D
  CONDITION eckit_EIGEN_FOUND
  SOURCES   test_Quad2D.cc
  LIBS      atlas
  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
)

ecbuild_add_test( TARGET atlas_test_Triag2D
  CONDITION eckit_EIGEN_FOUND
  SOURCES   test_Triag2D.cc
  LIBS      atlas
  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
)

ecbuild_add_test( TARGET atlas_test_interpolation_finite_element
  SOURCES   test_interpolation_finite_element.cc
  LIBS      atlas
  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
)

ecbuild_add_test( TARGET atlas_test_interpolation_unstructured_bilinear_lonlat
  SOURCES   test_interpolation_unstructured_bilinear_lonlat.cc
  LIBS      atlas
  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
)

ecbuild_add_test( TARGET atlas_test_interpolation_finite_element_cached
  SOURCES   test_interpolation_finite_element_cached.cc
  LIBS      atlas
  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
)

ecbuild_add_test( TARGET atlas_test_interpolation_grid_box_average
  SOURCES   test_interpolation_grid_box_average.cc
  LIBS      atlas
  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
)

ecbuild_add_test( TARGET atlas_test_interpolation_k_nearest_neighbours
  SOURCES   test_interpolation_k_nearest_neighbours.cc
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

ecbuild_add_test( TARGET atlas_test_interpolation_non_linear
  SOURCES  test_interpolation_non_linear.cc
  LIBS     atlas
  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
)

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

ecbuild_add_test( TARGET atlas_test_interpolation_structured2D_to_unstructured
  SOURCES  test_interpolation_structured2D_to_unstructured.cc
  LIBS     atlas
  MPI      2
  CONDITION eckit_HAVE_MPI
  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
)

ecbuild_add_test( TARGET atlas_test_interpolation_cubedsphere
  SOURCES  test_interpolation_cubedsphere.cc
  LIBS     atlas
  MPI      6
  CONDITION eckit_HAVE_MPI AND MPI_SLOTS GREATER_EQUAL 6 AND eckit_EIGEN_FOUND
  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
)
