# (C) Copyright 2013 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.


ecbuild_add_test( TARGET atlas_test_array_slicer
  SOURCES  test_array_slicer.cc
  LIBS     atlas
  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
)

ecbuild_add_test( TARGET atlas_test_array_foreach
  SOURCES  test_array_foreach.cc
  LIBS     atlas
  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
)

ecbuild_add_test( TARGET atlas_test_array
  SOURCES  test_array.cc
  LIBS     atlas
  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
)

#ecbuild_add_test( TARGET atlas_test_table
#  SOURCES  test_table.cc
#  LIBS     atlas
#  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
#)

ecbuild_add_test( TARGET atlas_test_array_vector
  SOURCES  test_array_vector.cc
  LIBS     atlas
  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
  CONDITION atlas_HAVE_GRIDTOOLS_STORAGE
)

ecbuild_add_test( TARGET atlas_test_svector
  SOURCES  test_svector.cc
  LIBS     atlas
  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}

)

ecbuild_add_test( TARGET atlas_test_array_view_util
  SOURCES  test_array_view_util.cc
  LIBS     atlas
  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
)

ecbuild_add_test( TARGET atlas_test_array_indexview
  SOURCES  test_indexview.cc
  LIBS     atlas
  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
)

if( CMAKE_BUILD_TYPE MATCHES "DEBUG" )
  set ( CMAKE_NVCC_FLAGS "-G" )
endif()

atlas_add_hic_test(
  TARGET atlas_test_array_kernel
  SOURCES test_array_kernel.hic
  LIBS atlas
  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
)
atlas_add_hic_test(
  TARGET atlas_test_vector_kernel
  SOURCES test_vector_kernel.hic
  LIBS atlas
  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
)

atlas_add_hic_test(
  TARGET atlas_test_svector_kernel
  SOURCES test_svector_kernel.hic
  LIBS atlas
  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
)

if(eckit_VERSION VERSION_GREATER_EQUAL 1.27.0)
# The test uses a header file from eckit "Overloaded.h" which only is available in v 1.27.0
ecbuild_add_test( TARGET atlas_test_array_view_variant
  SOURCES test_array_view_variant.cc
  LIBS atlas
  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
)
endif()
