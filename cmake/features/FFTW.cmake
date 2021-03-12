### FFTW ...

find_package(FFTW COMPONENTS double QUIET )
ecbuild_add_option( FEATURE FFTW
                    DESCRIPTION "Support for fftw"
                    REQUIRED_PACKAGES "FFTW COMPONENTS double QUIET" )
