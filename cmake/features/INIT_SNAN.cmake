### Init signaling NaN

if( ${CMAKE_BUILD_TYPE} MATCHES "Debug" )
  set( DEFAULT_INIT_SNAN ON )
else()
  set( DEFAULT_INIT_SNAN OFF )
endif()

ecbuild_add_option( FEATURE INIT_SNAN
                    DEFAULT ${DEFAULT_INIT_SNAN}
                    DESCRIPTION "Initialise atlas arrays with signaling_NaN (real types) or other invalid values (other types)" )

if( ${CMAKE_BUILD_TYPE} MATCHES "Debug" )
  if( NOT atlas_HAVE_INIT_SNAN )
    ecbuild_info( "Turning INIT_SNAN ON for Debug build" )
    set( atlas_HAVE_INIT_SNAN 1 )
  endif()
endif()

