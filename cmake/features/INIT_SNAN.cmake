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


cmake_push_check_state(RESET)
include(CheckSymbolExists)
  set(CMAKE_REQUIRED_DEFINITIONS -D_GNU_SOURCE)
  if(UNIX)
    set(CMAKE_REQUIRED_LIBRARIES m)
  endif()
  check_symbol_exists(feenableexcept "fenv.h" atlas_HAVE_FEENABLEEXCEPT)
  check_symbol_exists(fedisableexcept "fenv.h" atlas_HAVE_FEDISABLEEXCEPT)
  if( atlas_HAVE_FEENABLEEXCEPT )
      set( atlas_HAVE_FEENABLEEXCEPT 1 )
  else()
      set( atlas_HAVE_FEENABLEEXCEPT 0 )
  endif()
  if( atlas_HAVE_FEDISABLEEXCEPT )
      set( atlas_HAVE_FEDISABLEEXCEPT 1 )
  else()
      set( atlas_HAVE_FEDISABLEEXCEPT 0 )
  endif()
cmake_pop_check_state()
