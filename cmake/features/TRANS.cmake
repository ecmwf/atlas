### trans ...

set( atlas_HAVE_ECTRANS 0 )
if( ENABLE_TRANS OR NOT DEFINED ENABLE_TRANS )
    find_package( ectrans 1.1 COMPONENTS transi double QUIET )
    if( TARGET transi_dp )
        set( transi_FOUND TRUE )
        if( NOT TARGET transi )
            get_target_property( transi_dp_IMPORTED transi_dp IMPORTED )
            if( transi_dp_IMPORTED )
                set_target_properties( transi_dp PROPERTIES IMPORTED_GLOBAL TRUE) # required for aliasing imports
            endif()
            add_library( transi ALIAS transi_dp )
        endif()
        set( atlas_HAVE_ECTRANS 1 )
    else()
        find_package( transi 0.8 QUIET )
    endif()
endif()

ecbuild_add_option( FEATURE TRANS
                    DESCRIPTION "Support for IFS spectral transforms"
                    CONDITION transi_FOUND )
