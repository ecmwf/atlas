# (C) Copyright 2013 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

if( NOT DEFINED MPI_SLOTS )
  set( MPI_SLOTS 9999 )
endif()

######################################################
# Compiler bugs force us to disable some tests

if( CMAKE_Fortran_COMPILER_ID MATCHES "PGI"  AND
    CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 16.8 )
  set( atlas_fctest_parametrisation_DISABLED TRUE )
  ecbuild_warn( "PGI Fortran compiler version tested up to 16.7 has a compiler bug "
                "that leads to a segfaults for atlas_fctest_parametrisation. "
                "Disabling this test..." )
endif()

if( CMAKE_Fortran_COMPILER_ID MATCHES "Intel" )
  ## ATLAS-170 workaround internal compiler error

      ## Other workaround that makes ecbuild_remove_fortran_flags work when variables are empty
      string( TOUPPER ${CMAKE_BUILD_TYPE} CMAKE_BUILD_TYPE_CAPS )
      if( NOT CMAKE_Fortran_FLAGS )
          set( CMAKE_Fortran_FLAGS " " )
      endif()
      if( NOT CMAKE_Fortran_FLAGS_${CMAKE_BUILD_TYPE_CAPS} )
          set( CMAKE_Fortran_FLAGS_${CMAKE_BUILD_TYPE_CAPS} " " )
      endif()

  ecbuild_remove_fortran_flags( "-g" )
endif()

######################################################

# Macro atlas_add_cuda_test
# Envisioned to become part of ecbuild_add_test with a CUDA option
# in a future ecbuild release
macro( atlas_add_cuda_test )
  set( options           "" )
  set( single_value_args TARGET )
  set( multi_value_args  SOURCES LIBS ENVIRONMENT )
  cmake_parse_arguments( _PAR "${options}" "${single_value_args}" "${multi_value_args}"  ${_FIRST_ARG} ${ARGN} )

  if(_PAR_UNPARSED_ARGUMENTS)
    ecbuild_critical("Unknown keywords given to atlas_add_cuda_test(): \"${_PAR_UNPARSED_ARGUMENTS}\"")
  endif()

  if( HAVE_CUDA AND HAVE_TESTS )

    ecbuild_debug("atlas_add_cuda_test: Adding test ${_PAR_TARGET}")

    list( APPEND _PAR_ENVIRONMENT "ATLAS_RUN_NGPUS=1" )
    list( REMOVE_DUPLICATES _PAR_ENVIRONMENT )

    list( APPEND _libs ${_PAR_LIBS} )
    if( _libs )
      ecbuild_debug("atlas_add_cuda_test: Test ${_PAR_TARGET} explicitly links with libraries ${_libs}")
    endif()

    if( atlas_CUDA_LANGUAGE_ENABLED )

      ecbuild_add_test( TARGET      ${_PAR_TARGET}
                        SOURCES     ${_PAR_SOURCES}
                        LIBS        ${_libs}
                        ENVIRONMENT ${_PAR_ENVIRONMENT} )

    else()

      cuda_add_executable (${_PAR_TARGET} ${_PAR_SOURCES})

      if( _libs )
        target_link_libraries(${_PAR_TARGET} ${_libs})
      endif()

      # whatever project settings are, we always build tests with the build_rpath, not the install_rpath
      set_property( TARGET ${_PAR_TARGET} PROPERTY BUILD_WITH_INSTALL_RPATH FALSE )
      set_property( TARGET ${_PAR_TARGET} PROPERTY SKIP_BUILD_RPATH         FALSE )

      if( CMAKE_CROSSCOMPILING_EMULATOR )
        add_test ( NAME ${_PAR_TARGET} COMMAND ${CMAKE_CROSSCOMPILING_EMULATOR} ${_PAR_TARGET} )
      else()
        add_test ( NAME ${_PAR_TARGET} COMMAND ${_PAR_TARGET} )
      endif()

      set_property( TEST ${_PAR_TARGET} APPEND PROPERTY ENVIRONMENT "${_PAR_ENVIRONMENT}" )

    endif()
  endif()
endmacro()

if( ECKIT_INCLUDE_DIRS ) # eckit not yet ported to CMake3
  include_directories( ${ECKIT_INCLUDE_DIRS} )
endif()

if( NOT DEFINED transi_VERSION )
  if( TRANSI_VERSION VERSION_LESS 0.6 )# transi not yet ported to CMake3
    include_directories( ${TRANSI_INCLUDE_DIRS} )
  endif()
endif()

if( HAVE_CUDA )
  set( ATLAS_TEST_ENVIRONMENT "ATLAS_RUN_NGPUS=1" )
endif()

add_subdirectory( util )
add_subdirectory( runtime )

if (atlas_HAVE_ATLAS_GRID)
  add_subdirectory( projection )
  add_subdirectory( grid )
endif()

if (atlas_HAVE_ATLAS_FIELD)
  add_subdirectory( acc )
  add_subdirectory( array )
  add_subdirectory( field )
endif()

if (atlas_HAVE_ATLAS_FUNCTIONSPACE)
  add_subdirectory( parallel )
  add_subdirectory( mesh )
  add_subdirectory( functionspace )
  add_subdirectory( io )
  add_subdirectory( output )
  add_subdirectory( numerics )
  add_subdirectory( trans )
  add_subdirectory( interpolation )
  add_subdirectory( linalg )
  add_subdirectory( redistribution )
endif()
add_subdirectory( acceptance_tests )
add_subdirectory( export_tests )
