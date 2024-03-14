if( atlas_HAVE_ATLAS_FIELD AND NOT HAVE_GRIDTOOLS_STORAGE )

### Native_Storage storage module

    ### Native_Storage may search for CUDA, which searches for "Threads"
    ### Set THREADS_HAVE_PTHREAD_ARG variable to false so that it can be recomputed based on
    ### THREADS_PREFER_PTHREAD_FLAG, in case other project had it on a different setting.
    ### This is certainly a CMake bug ( see ECKIT-426 )
    set( THREADS_HAVE_PTHREAD_ARG FALSE )
    if( NOT DEFINED THREADS_PREFER_PTHREAD_FLAG )
      set( THREADS_PREFER_PTHREAD_FLAG 1 )
    endif()

    #ecbuild_add_option(
    # FEATURE NATIVE_STORAGE
    #DESCRIPTION "Arrays internally use Native storage layer"
    #)

ecbuild_add_option( FEATURE CUDA
                    DESCRIPTION "Enable CUDA support via Native CUDA backend"
                  )
set( ATLAS_NATIVE_STORAGE_BACKEND_HOST 0 )
set( ATLAS_NATIVE_STORAGE_BACKEND_CUDA 0 )
set( atlas_HAVE_NATIVE_STORAGE 1 ) # always true

if( atlas_HAVE_CUDA )

    ecbuild_info( "Native storage with CUDA support" )

    # Logic to check if we can use enable_language( CUDA )
    #   - CMake already supports it (as Native_Storage requires version that supports it)
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

    find_package( CUDAToolkit QUIET)

    #    target_include_directories( atlas INTERFACE ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES} )
    #target_link_libraries( atlas INTERFACE ${CUDA_CUDART_LIBRARY} )

    set( ATLAS_NATIVE_STORAGE_BACKEND_CUDA 1 )

else()

    set( ATLAS_NATIVE_STORAGE_BACKEND_HOST 1 )

endif()

else()
    set( ATLAS_NATIVE_STORAGE_BACKEND_CUDA 0 )
    set( ATLAS_NATIVE_STORAGE_BACKEND_HOST 1 )
endif()
