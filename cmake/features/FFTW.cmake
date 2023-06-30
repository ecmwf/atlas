### FFTW ...

if( atlas_HAVE_ATLAS_TRANS )

ecbuild_add_option( FEATURE FFTW
                    DESCRIPTION "Support for fftw"
                    REQUIRED_PACKAGES "FFTW COMPONENTS double QUIET" )

if( NOT HAVE_FFTW )
    unset( FFTW_LIBRARIES )
    unset( FFTW_INCLUDES )
endif()

endif()