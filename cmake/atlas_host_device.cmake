# (C) Copyright 2013 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

function( create_cuda_wrapper variable )
  set( options "" )
  set( single_value_args SOURCE )
  set( multi_value_args "" )
  cmake_parse_arguments( _PAR "${options}" "${single_value_args}" "${multi_value_args}"  ${_FIRST_ARG} ${ARGN} )

  get_filename_component(directory ${_PAR_SOURCE} DIRECTORY)
  get_filename_component(base      ${_PAR_SOURCE} NAME_WE)
  get_filename_component(name      ${_PAR_SOURCE} NAME)
  get_filename_component(abspath   ${_PAR_SOURCE} ABSOLUTE)

  if( directory )
    set(cuda_wrapper ${CMAKE_CURRENT_BINARY_DIR}/${directory}/${base}.cu)
  else()
    set(cuda_wrapper ${CMAKE_CURRENT_BINARY_DIR}/${base}.cu)
  endif()
  set(${variable} ${cuda_wrapper} PARENT_SCOPE)
  set(content
"
#include \"atlas/${directory}/${name}\"
")
  if( ${abspath} IS_NEWER_THAN ${cuda_wrapper} )
    file(WRITE ${cuda_wrapper} "${content}")
  endif()
endfunction()



function( atlas_host_device srclist )
  set( options "" )
  set( single_value_args "" )
  set( multi_value_args SOURCES )
  cmake_parse_arguments( _PAR "${options}" "${single_value_args}" "${multi_value_args}" ${ARGN} )

  if( ATLAS_GRIDTOOLS_STORAGE_BACKEND_CUDA )
    set( use_cuda_srclist ${_PAR_SOURCES} )

    foreach( src ${use_cuda_srclist} )
      create_cuda_wrapper( cuda_wrapper SOURCE ${src} )
      list( APPEND ${srclist} ${cuda_wrapper} )
      ecbuild_list_exclude_pattern( LIST ${srclist} REGEX ${src} )
    endforeach()

    set( ${srclist} "${${srclist}}" PARENT_SCOPE )
  endif()
endfunction()
