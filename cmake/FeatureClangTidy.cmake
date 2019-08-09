set( CLANG_TIDY_SUPPORTED_COMPILERS     GNU Clang )

set( CLANG_TIDY_DEFAULT OFF )
foreach( _clang_tidy_supported_compiler ${CLANG_TIDY_SUPPORTED_COMPILERS} )
  if( CMake_CXX_COMPILER_ID MATCHES _clang_tidy_supported_compiler )
    set( CLANG_TIDY_DEFAULT ON )
  endif()
endforeach()

find_program( CLANG_TIDY_EXE NAMES "clang-tidy" )
if( CLANG_TIDY_EXE )
    ecbuild_info( "Found clang-tidy: ${CLANG_TIDY_EXE}" )
    if( NOT CLANG_TIDY_DEFAULT AND NOT DEFINED ENABLE_CLANG_TIDY )
      ecbuild_info( "Feature CLANG_TIDY default OFF for compiler ${CMAKE_CXX_COMPILER_ID}" )
    endif()
endif()

ecbuild_add_option( FEATURE CLANG_TIDY
                    DEFAULT ${CLANG_TIDY_DEFAULT}
                    DESCRIPTION "Use clang-tidy"
                    CONDITION CLANG_TIDY_EXE )

if (HAVE_CLANG_TIDY)
  set(CLANG_TIDY_CHECKS "-*,readability-braces-around-statements,redundant-string-init")
  set(CMAKE_CXX_CLANG_TIDY "${CLANG_TIDY_EXE};-checks=${CLANG_TIDY_CHECKS};-header-filter='${CMAKE_SOURCE_DIR}/*'")
endif()
