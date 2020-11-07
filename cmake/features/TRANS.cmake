### trans ...

if( ENABLE_TRANS OR NOT DEFINED ENABLE_TRANS )
    find_package( trans 47.2 COMPONENTS transi double QUIET )
    if( TARGET transi_dp )
        set( transi_FOUND TRUE )
        if( NOT TARGET transi )
            set_target_properties( transi_dp PROPERTIES IMPORTED_GLOBAL TRUE) # required for aliasing imports
            add_library( transi ALIAS transi_dp )
        endif()
    else()
        find_package( transi 0.8 QUIET )
    endif()
endif()

ecbuild_add_option( FEATURE TRANS
                    DESCRIPTION "Support for IFS spectral transforms"
                    CONDITION transi_FOUND )
