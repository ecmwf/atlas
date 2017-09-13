# (C) Copyright 1996-2017 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

ecbuild_add_test( TARGET atlas_test_array
  SOURCES  test_array.cc
  LIBS     atlas
)

#ecbuild_add_test( TARGET atlas_test_table
#  SOURCES  test_table.cc
#  LIBS     atlas
#)

if( CMAKE_BUILD_TYPE MATCHES "DEBUG" )
  set ( CMAKE_NVCC_FLAGS "-G" )
endif()

atlas_add_cuda_test(
  TARGET atlas_test_array_kernel
  BOOST
  SOURCES test_array_kernel.cu
  LIBS atlas
)
atlas_add_cuda_test(
  TARGET atlas_test_vector_kernel
  BOOST
  SOURCES test_vector_kernel.cu
  LIBS atlas
)

