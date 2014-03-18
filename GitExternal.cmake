# configures an external git repository
# Usage:
# git_external(<directory> <giturl> <gittag> [RESET <files>])
# git_external_manage(<file>)
# is used for a managed .gitexternal file which is parsed and
# updated by this function, which calls git_external and provides
# an update target to bump the tag to the master revision by
# recreating the given file. The file has to contain:
# include(GitExternal)
# git_external_manage(${CMAKE_CURRENT_LIST_FILE})
# #-> CMake/common https://github.com/Eyescale/CMake.git 8324611
# Where the last line can appear multiple times, once for each
# external repo It has to start with '#->'. This function parses
# the file for this pattern and then calls git_external on each
# found entry. The entries are encoded into a target which will
# recreate the file on invocation. Since only the '#-> ' lines are
# parsed, everything else is ignored. Since the file is under git
# control, there should be no risk of data loss.

find_package(Git)
if(NOT GIT_EXECUTABLE)
  return()
endif()

include(CMakeParseArguments)

function(GIT_EXTERNAL_MANAGE FILE)
  if(NOT EXISTS ${FILE})
    message(FATAL_ERROR "Can't open ${FILE}")
  endif()

  file(READ ${FILE} GIT_EXTERNAL_FILE)
  string(REGEX REPLACE "\n" ";" GIT_EXTERNAL_FILE "${GIT_EXTERNAL_FILE}")
  foreach(LINE ${GIT_EXTERNAL_FILE})
    string(REGEX REPLACE "^#->[ ]*(.+[ ]+.+[ ]+.+)$" "\\1" DATA ${LINE})
    if(NOT LINE STREQUAL DATA)
      string(REGEX REPLACE "[ ]+" ";" DATA "${DATA}")
      list(LENGTH DATA DATA_LENGTH)
      if(DATA_LENGTH EQUAL 3)
        list(GET DATA 0 DIR)
        list(GET DATA 1 REPO)
        list(GET DATA 2 TAG)

        # pull in identified external
        git_external(${DIR} ${REPO} ${TAG})

        # Create update script and target to bump external spec
        if(NOT TARGET update)
          add_custom_target(update)
        endif()
        if(NOT TARGET update_git_external)
          add_custom_target(update_git_external)
          add_dependencies(update update_git_external)
        endif()
        if(NOT TARGET update_git_external_header)
          set(GIT_EXTERNAL_SCRIPT
            "${CMAKE_CURRENT_BINARY_DIR}/gitupdateexternal.cmake")
          file(WRITE "${GIT_EXTERNAL_SCRIPT}"
            "file(WRITE ${FILE} \"# -*- mode: cmake -*-
include(GitExternal)
git_external_manage(\\\${CMAKE_CURRENT_LIST_FILE})
\")")
          add_custom_target(update_git_external_header
            COMMAND ${CMAKE_COMMAND} -P ${GIT_EXTERNAL_SCRIPT}
            COMMENT "Recreate ${FILE}"
            WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")
        endif()

        string(REPLACE "/" "_" GIT_EXTERNAL_NAME ${DIR}) # unique, flat name
        set(GIT_EXTERNAL_SCRIPT
          "${CMAKE_CURRENT_BINARY_DIR}/gitupdate${GIT_EXTERNAL_NAME}.cmake")
        file(WRITE "${GIT_EXTERNAL_SCRIPT}" "
execute_process(COMMAND ${GIT_EXECUTABLE} fetch --all -q
WORKING_DIRECTORY ${DIR})
execute_process(
COMMAND ${GIT_EXECUTABLE} show-ref --hash=7 refs/remotes/origin/master
OUTPUT_VARIABLE newref WORKING_DIRECTORY ${DIR})
if(newref)
file(APPEND ${FILE} \"#-> ${DIR} ${REPO} \${newref}\n\")
else()
  file(APPEND ${FILE} \"#-> ${DIR} ${REPO} ${TAG}\n\")
endif()")
        add_custom_target(update_git_external_${GIT_EXTERNAL_NAME}
          COMMAND ${CMAKE_COMMAND} -P ${GIT_EXTERNAL_SCRIPT}
          COMMENT "Update ${REPO} in ${FILE}"
          DEPENDS update_git_external_header
          WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")
        add_dependencies(update_git_external
          update_git_external_${GIT_EXTERNAL_NAME})
      endif()
    endif()
  endforeach()
endfunction()

include(CMakeParseArguments)

macro( GIT_EXTERNAL )

    set( single_value_args DIR REPO TAG HISTORY )

    cmake_parse_arguments( _PAR "${options}" "${single_value_args}" "${multi_value_args}" ${_FIRST_ARG} ${ARGN} )

    if(_PAR_UNPARSED_ARGUMENTS)
      message(FATAL_ERROR "Unknown keywords given to git_external(): \"${_PAR_UNPARSED_ARGUMENTS}\"")
    endif()

    # default is no history
    if( NOT DEFINED _PAR_HISTORY )
      set( _PAR_HISTORY 0 )
    endif()

    if( NOT ${_PAR_HISTORY} )
      set( depth "--depth=1")
    endif()

#function(GIT_EXTERNAL DIR REPO TAG)
#  cmake_parse_arguments(GIT_EXTERNAL "" "" "RESET" ${ARGN})

  get_filename_component(_PAR_DIR "${_PAR_DIR}" ABSOLUTE)
  get_filename_component(GIT_EXTERNAL_DIR "${_PAR_DIR}/.." ABSOLUTE)

  if(NOT EXISTS "${_PAR_DIR}")
    message(STATUS "git clone ${_PAR_REPO} ${_PAR_DIR} ${depth}")
    execute_process(
      COMMAND "${GIT_EXECUTABLE}" clone "${_PAR_REPO}" "${_PAR_DIR}" "${depth}"
      RESULT_VARIABLE nok ERROR_VARIABLE error
      WORKING_DIRECTORY "${GIT_EXTERNAL_DIR}")
    if(nok)
      message(FATAL_ERROR "${_PAR_DIR} git clone failed: ${error}\n")
    endif()
  endif()

  if(IS_DIRECTORY "${_PAR_DIR}/.git")
    foreach(GIT_EXTERNAL_RESET_FILE ${GIT_EXTERNAL_RESET})
      execute_process(
        COMMAND "${GIT_EXECUTABLE}" reset -q "${GIT_EXTERNAL_RESET_FILE}"
        RESULT_VARIABLE nok ERROR_VARIABLE error
        WORKING_DIRECTORY "${_PAR_DIR}")
      execute_process(
        COMMAND "${GIT_EXECUTABLE}" checkout -q -- "${GIT_EXTERNAL_RESET_FILE}"
        RESULT_VARIABLE nok ERROR_VARIABLE error
        WORKING_DIRECTORY "${_PAR_DIR}")
    endforeach()

    execute_process(COMMAND "${GIT_EXECUTABLE}" fetch --all -q
      RESULT_VARIABLE nok ERROR_VARIABLE error
      WORKING_DIRECTORY "${_PAR_DIR}")
    if(nok)
      message(STATUS "Update of ${_PAR_DIR} failed:\n ${error}")
    endif()

    execute_process(
      COMMAND "${GIT_EXECUTABLE}" checkout -q "${_PAR_TAG}"
      RESULT_VARIABLE nok ERROR_VARIABLE error
      WORKING_DIRECTORY "${_PAR_DIR}"
      )
    if(nok)
      message(FATAL_ERROR "${_PAR_DIR} git checkout ${_PAR_TAG} failed: ${error}\n")
    endif()
  endif()
#endfunction()
endmacro()

