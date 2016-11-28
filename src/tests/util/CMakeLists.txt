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
  CONDITION ATLAS_HAVE_GRIDTOOLS_STORAGE
)

ecbuild_add_test( TARGET atlas_test_metadata
  BOOST
  MPI      4
  SOURCES  test_metadata.cc
  LIBS     atlas
)

# Macro ecbuild_add_cuda_test
# Envisioned to become part of ecbuild_add_test with a CUDA option
# in a future ecbuild release
macro( ecbuild_add_cuda_test )
  set( options           BOOST )
  set( single_value_args TARGET )
  set( multi_value_args  SOURCES LIBS )
  cmake_parse_arguments( _PAR "${options}" "${single_value_args}" "${multi_value_args}"  ${_FIRST_ARG} ${ARGN} )

  if(_PAR_UNPARSED_ARGUMENTS)
    ecbuild_critical("Unknown keywords given to add_cuda_test(): \"${_PAR_UNPARSED_ARGUMENTS}\"")
  endif()

  if( ENABLE_GPU AND ENABLE_TESTS )

    list( APPEND _libs ${_PAR_LIBS} )

    if( _PAR_BOOST )
      list( APPEND _libs ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY} ${Boost_TEST_EXEC_MONITOR_LIBRARY} )
    endif()

    cuda_add_executable (${_PAR_TARGET} ${_PAR_SOURCES})
    if( _libs )
      target_link_libraries(${_PAR_TARGET} ${_libs})
    endif()

    # whatever project settings are, we always build tests with the build_rpath, not the install_rpath
    set_property( TARGET ${_PAR_TARGET} PROPERTY BUILD_WITH_INSTALL_RPATH FALSE )
    set_property( TARGET ${_PAR_TARGET} PROPERTY SKIP_BUILD_RPATH         FALSE )

    add_test (${_PAR_TARGET} ${_PAR_TARGET})

  endif()
endmacro()

#TODO Do this only in DEBUG mode
#set ( CMAKE_NVCC_FLAGS "-G" )
ecbuild_add_cuda_test(
  TARGET atlas_test_array_kernel
  BOOST
  SOURCES test_array_kernel.cu
  LIBS atlas
)
