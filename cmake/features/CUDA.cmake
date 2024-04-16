
ecbuild_add_option( FEATURE CUDA
                    DESCRIPTION "Enable CUDA support"
                    DEFAULT OFF
                  )

if( HAVE_CUDA )

    enable_language( CUDA )
    ecbuild_info( "CUDA language enabled" )

    find_package( CUDAToolkit REQUIRED )

endif()

