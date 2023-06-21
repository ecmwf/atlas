### FFTW ...

ecbuild_add_option( FEATURE FFTW
                    DESCRIPTION "Support for fftw"
                    CONDITION atlas_HAVE_ATLAS_NUMERICS
                    REQUIRED_PACKAGES "FFTW COMPONENTS double QUIET" )

if( NOT HAVE_FFTW )
    unset( FFTW_LIBRARIES )
    unset( FFTW_INCLUDES )
endif()
