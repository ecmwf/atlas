# (C) Copyright 2019 ECMWF.
#
# This file is covered by the LICENSING file in the root of this project.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

function( run_command COMMAND )
  set( COMMAND ${ARGV} )
  string(REPLACE ";" " " PRINT_COMMAND "${COMMAND}" )
  message( "cd ${CMAKE_CURRENT_BINARY_DIR} && ${PRINT_COMMAND}" )
  execute_process(
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
      COMMAND ${COMMAND}
      RESULT_VARIABLE res
  )
  if(res) 
    message(FATAL_ERROR "[ ${COMMAND} ] failed")
  endif() 
endfunction()

message( "Working directory: ${CMAKE_CURRENT_BINARY_DIR}")

list( GET FILES 0 ref )
foreach( file ${FILES} )
  if( EXISTS ${file} )
    if( NOT ${file} STREQUAL ${ref} )
      run_command( ${CMAKE_COMMAND} -E compare_files ${ref} ${file} )
    endif()
  else()
    message( STATUS "Skipping ${file} as it does not exist" )
  endif()
endforeach()
