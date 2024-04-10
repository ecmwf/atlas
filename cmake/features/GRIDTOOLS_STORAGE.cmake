if( atlas_HAVE_ATLAS_FIELD AND (ENABLE_GRIDTOOLS_STORAGE OR atlas_ENABLE_GRIDTOOLS_STORAGE) )

### GridTools storage module

    ### GridTools may search for CUDA, which searches for "Threads"
    ### Set THREADS_HAVE_PTHREAD_ARG variable to false so that it can be recomputed based on
    ### THREADS_PREFER_PTHREAD_FLAG, in case other project had it on a different setting.
    ### This is certainly a CMake bug ( see ECKIT-426 )
    set( THREADS_HAVE_PTHREAD_ARG FALSE )
    if( NOT DEFINED THREADS_PREFER_PTHREAD_FLAG )
      set( THREADS_PREFER_PTHREAD_FLAG 1 )
    endif()

find_package( GridTools QUIET 
    HINTS ${GridTools_ROOT}/lib/cmake
          $ENV{GridTools_ROOT}/lib/cmake
          ${CMAKE_PREFIX_PATH}/lib/cmake
          ${CMAKE_INSTALL_PREFIX}/lib/cmake
          ${GridTools_BINARY_DIR} )
ecbuild_add_option(
  FEATURE GRIDTOOLS_STORAGE
  DESCRIPTION "Arrays internally use GridTools storage layer"
  CONDITION GridTools_FOUND )

set( ATLAS_GRIDTOOLS_STORAGE_BACKEND_HOST 0 )
set( ATLAS_GRIDTOOLS_STORAGE_BACKEND_CUDA 0 )

if( atlas_HAVE_GRIDTOOLS_STORAGE )
  if( GRIDTOOLS_HAS_BACKEND_CUDA )
    if( atlas_HAVE_CUDA )
      ecbuild_info( "GridTools found with CUDA support -> backend CUDA" )
      set( ATLAS_GRIDTOOLS_STORAGE_BACKEND_CUDA 1 )
    else()
      ecbuild_info( "GridTools found with CUDA support, but atlas does not have CUDA enabled -> backend HOST" )
      set( ATLAS_GRIDTOOLS_STORAGE_BACKEND_HOST 1 )
    endif()
  else()
    ecbuild_info( "GridTools found without CUDA support -> backend HOST" )
    set( ATLAS_GRIDTOOLS_STORAGE_BACKEND_HOST 1 )
  endif()
endif()

else()
  set( HAVE_GRIDTOOLS_STORAGE 0 )
  set( atlas_HAVE_GRIDTOOLS_STORAGE 0 )
  set( ATLAS_GRIDTOOLS_STORAGE_BACKEND_HOST 0 )
  set( ATLAS_GRIDTOOLS_STORAGE_BACKEND_CUDA 0 )
endif()

