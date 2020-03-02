### Bounds checking
if( ${CMAKE_BUILD_TYPE} MATCHES "Debug" )
  set( DEFAULT_BOUNDSCHECKING ON )
else()
  set( DEFAULT_BOUNDSCHECKING OFF )
endif()

ecbuild_add_option( FEATURE BOUNDSCHECKING
                    DEFAULT ${DEFAULT_BOUNDSCHECKING}
                    DESCRIPTION "Bounds checking for atlas::ArrayView and atlas::IndexView" )


if( ${CMAKE_BUILD_TYPE} MATCHES "Debug" )
  if( NOT atlas_HAVE_BOUNDSCHECKING )
    ecbuild_info( "Turning BOUNDSCHECKING ON for Debug build" )
    set( atlas_HAVE_BOUNDSCHECKING 1 )
  endif()
endif()

