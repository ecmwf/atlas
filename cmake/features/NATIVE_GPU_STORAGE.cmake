if( atlas_HAVE_ATLAS_FIELD AND ENABLE_NATIVE_GPU_STORAGE )

### Native_GPU_Storage storage module

    ### Native_GPU_Storage may search for CUDA, which searches for "Threads"
    ### Set THREADS_HAVE_PTHREAD_ARG variable to false so that it can be recomputed based on
    ### THREADS_PREFER_PTHREAD_FLAG, in case other project had it on a different setting.
    ### This is certainly a CMake bug ( see ECKIT-426 )
    set( THREADS_HAVE_PTHREAD_ARG FALSE )
    if( NOT DEFINED THREADS_PREFER_PTHREAD_FLAG )
      set( THREADS_PREFER_PTHREAD_FLAG 1 )
    endif()

ecbuild_add_option(
  FEATURE NATIVE_GPU_STORAGE
  DESCRIPTION "Arrays internally use Native_GPU_Storage storage layer"
)

ecbuild_add_option( FEATURE CUDA
                    DESCRIPTION "Enable CUDA support via Native_GPU_Storage CUDA backend"
                    )

set( ATLAS_NATIVE_GPU_STORAGE_BACKEND_HOST 0 )
set( ATLAS_NATIVE_GPU_STORAGE_BACKEND_CUDA 0 )
set( atlas_HAVE_NATIVE_GPU_STORAGE 1 )

if( atlas_HAVE_NATIVE_GPU_STORAGE )
    if( atlas_HAVE_CUDA )

        ecbuild_info( "Native GPU Storage with CUDA support" )

        # Logic to check if we can use enable_language( CUDA )
        #   - CMake already supports it (as Native_GPU_Storage requires version that supports it)
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

        set( ATLAS_NATIVE_GPU_STORAGE_BACKEND_CUDA 1 )

    else()

        set( ATLAS_NATIVE_GPU_STORAGE_BACKEND_HOST 1 )

    endif()
endif()

else()
    set( HAVE_NATIVE_GPU_STORAGE 1 )
    set( atlas_HAVE_NATIVE_GPU_STORAGE 0 )
    set( ATLAS_NATIVE_GPU_STORAGE_BACKEND_CUDA 0 )
    set( ATLAS_NATIVE_GPU_STORAGE_BACKEND_HOST 0 )
endif()
