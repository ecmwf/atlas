### trans ...

ecbuild_find_package( NAME transi VERSION 0.4.4 QUIET )
ecbuild_add_option( FEATURE TRANS
                    DESCRIPTION "Support for spectral transforms"
                    CONDITION transi_FOUND )
