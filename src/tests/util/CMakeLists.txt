# (C) Copyright 1996-2016 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

if( HAVE_FCTEST )

  add_fctest( TARGET atlas_fctest_refcounting
    LINKER_LANGUAGE Fortran
    SOURCES         fctest_refcounting.F90
    LIBS            atlas_f
  )

  add_fctest( TARGET atlas_fctest_logging
    MPI             1
    LINKER_LANGUAGE Fortran
    SOURCES         fctest_logging.F90
    LIBS            atlas_f
  )

  add_fctest( TARGET atlas_fctest_error
    MPI             1
    LINKER_LANGUAGE Fortran
    SOURCES         fctest_error.F90
    LIBS            atlas_f
  )

  add_fctest( TARGET atlas_fctest_value
    MPI             1
    LINKER_LANGUAGE Fortran
    SOURCES         fctest_value.F90
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

ecbuild_add_test( TARGET atlas_test_array
  BOOST
  SOURCES  test_array.cc
  LIBS     atlas
)

ecbuild_add_test( TARGET atlas_test_metadata
  BOOST
  MPI      4
  SOURCES  test_metadata.cc
  LIBS     atlas
)

if( ENABLE_GPU )
ADD_LIBRARY(boost_unit_test_framework STATIC IMPORTED)
cuda_add_executable (test_array_kernel test_array_kernel.cu)
  target_link_libraries(test_array_kernel eckit ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY} cudart pthread)
#  set_target_properties(test_array_kernel PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${DIR} OUTPUT_NAME ${unit_test})
  add_test (test_array_kernel test_array_kernel)
endif()

