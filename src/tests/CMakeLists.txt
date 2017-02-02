# (C) Copyright 1996-2017 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.



# Macro atlas_add_cuda_test
# Envisioned to become part of ecbuild_add_test with a CUDA option
# in a future ecbuild release
macro( atlas_add_cuda_test )
  set( options           BOOST )
  set( single_value_args TARGET )
  set( multi_value_args  SOURCES LIBS )
  cmake_parse_arguments( _PAR "${options}" "${single_value_args}" "${multi_value_args}"  ${_FIRST_ARG} ${ARGN} )

  if(_PAR_UNPARSED_ARGUMENTS)
    ecbuild_critical("Unknown keywords given to atlas_add_cuda_test(): \"${_PAR_UNPARSED_ARGUMENTS}\"")
  endif()

  if( ATLAS_HAVE_GRIDTOOLS_STORAGE AND ENABLE_GPU AND ENABLE_TESTS )

    ecbuild_debug("atlas_add_cuda_test: Adding test ${_PAR_TARGET}")

    list( APPEND _libs ${_PAR_LIBS} )

    if( _PAR_BOOST )
      list( APPEND _libs ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY} ${Boost_TEST_EXEC_MONITOR_LIBRARY} )
    endif()

    cuda_add_executable (${_PAR_TARGET} ${_PAR_SOURCES})
    if( _libs )
      ecbuild_debug("atlas_add_cuda_test: Test ${_PAR_TARGET} explicitly links with libraries ${_libs}")
      target_link_libraries(${_PAR_TARGET} ${_libs})
    endif()

    # whatever project settings are, we always build tests with the build_rpath, not the install_rpath
    set_property( TARGET ${_PAR_TARGET} PROPERTY BUILD_WITH_INSTALL_RPATH FALSE )
    set_property( TARGET ${_PAR_TARGET} PROPERTY SKIP_BUILD_RPATH         FALSE )

    add_test (${_PAR_TARGET} ${_PAR_TARGET})

  endif()
endmacro()


add_subdirectory( array )
add_subdirectory( util )
add_subdirectory( parallel )
add_subdirectory( field )
add_subdirectory( grid )
add_subdirectory( mesh )
add_subdirectory( functionspace )
add_subdirectory( io )
add_subdirectory( numerics )
add_subdirectory( trans )
add_subdirectory( interpolation )

