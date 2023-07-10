if( atlas_HAVE_ATLAS_FIELD )

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

ecbuild_add_option( FEATURE CUDA
                    DESCRIPTION "Enable CUDA support via GridTools CUDA backend"
                    CONDITION GRIDTOOLS_HAS_BACKEND_CUDA )

set( ATLAS_GRIDTOOLS_STORAGE_BACKEND_HOST 0 )
set( ATLAS_GRIDTOOLS_STORAGE_BACKEND_CUDA 0 )

if( atlas_HAVE_GRIDTOOLS_STORAGE )

  if( atlas_HAVE_CUDA )

    ecbuild_info( "GridTools found with CUDA support" )

    # Logic to check if we can use enable_language( CUDA )
    #   - CMake already supports it (as GridTools requires version that supports it)
    #   - ecbuild version 3.3 added support
    #   - overriding mechanism possible with cached "atlas_CUDA_LANGUAGE_ENABLED" variable
    if( DEFINED ecbuild_VERSION AND NOT ecbuild_VERSION VERSION_LESS 3.3 )
      set( atlas_CUDA_LANGUAGE_ENABLED_DEFAULT ON )
    else()
      set( atlas_CUDA_LANGUAGE_ENABLED_DEFAULT OFF )
    endif()
    set( atlas_CUDA_LANGUAGE_ENABLED ${atlas_CUDA_LANGUAGE_ENABLED_DEFAULT} CACHE BOOL "atlas enables CUDA language" )

    if( atlas_CUDA_LANGUAGE_ENABLED )
      enable_language( CUDA )
      ecbuild_info( "CUDA language enabled" )
    else()
      ecbuild_info("CUDA enabled through find_package(CUDA) instead of enable_language(CUDA)")
      find_package( CUDA )
    endif()

    set( ATLAS_GRIDTOOLS_STORAGE_BACKEND_CUDA 1 )

  else()

    set( ATLAS_GRIDTOOLS_STORAGE_BACKEND_HOST 1 )

  endif()

endif()

else()
  set( HAVE_GRIDTOOLS_STORAGE 1 )
  set( atlas_HAVE_GRIDTOOLS_STORAGE 0 )
  set( ATLAS_GRIDTOOLS_STORAGE_BACKEND_CUDA 0 )
  set( ATLAS_GRIDTOOLS_STORAGE_BACKEND_HOST 0 )
endif()