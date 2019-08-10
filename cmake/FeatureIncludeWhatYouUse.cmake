find_program( INCLUDE_WHAT_YOU_USE_EXE NAMES "include-what-you-use" )
if( INCLUDE_WHAT_YOU_USE_EXE )
    ecbuild_info( "Found include-what-you-use: ${INCLUDE_WHAT_YOU_USE_EXE}" )
endif()

ecbuild_add_option( FEATURE INCLUDE_WHAT_YOU_USE
                    DEFAULT OFF # Need clang compiler?
                    DESCRIPTION "Use include-what-you-use clang-tool"
                    CONDITION INCLUDE_WHAT_YOU_USE_EXE )

if( HAVE_INCLUDE_WHAT_YOU_USE )
    set( CMAKE_CXX_INCLUDE_WHAT_YOU_USE ${INCLUDE_WHAT_YOU_USE_EXE} )
endif()
