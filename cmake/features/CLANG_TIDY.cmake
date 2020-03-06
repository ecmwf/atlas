set( CLANG_TIDY_SUPPORTED_COMPILERS     GNU Clang )

set( CLANG_TIDY_DEFAULT OFF )
foreach( _clang_tidy_supported_compiler ${CLANG_TIDY_SUPPORTED_COMPILERS} )
  if( CMAKE_CXX_COMPILER_ID MATCHES ${_clang_tidy_supported_compiler} )
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

  # Uncomment to apply fixes. Make sure to use a clean build, and apply clang-format afterwards!
  # set( CLANG_TIDY_FIXIT ";-fix" )

  set( CLANG_TIDY_CHECKS "-*" )
  foreach( _clang_tidy_check
    readability-braces-around-statements
    redundant-string-init
    modernize-use-nullptr
    modernize-use-using
    modernize-use-override
    modernize-use-emplace
    modernize-use-equals-default
    modernize-use-equals-delete
    )
    set( CLANG_TIDY_CHECKS "${CLANG_TIDY_CHECKS},${_clang_tidy_check}" )
  endforeach()
  set( CMAKE_CXX_CLANG_TIDY "${CLANG_TIDY_EXE};-checks=${CLANG_TIDY_CHECKS};-header-filter='${CMAKE_SOURCE_DIR}/*'${CLANG_TIDY_FIXIT}" )
endif()
