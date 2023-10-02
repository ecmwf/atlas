if( atlas_HAVE_ATLAS_TRANS )

### trans ...

if( NOT DEFINED ATLAS_ENABLE_ECTRANS AND DEFINED ATLAS_ENABLE_TRANS )
  ecbuild_warn("Atlas option ATLAS_ENABLE_TRANS is deprecated in favour of ATLAS_ENABLE_ECTRANS")
  set( ATLAS_ENABLE_ECTRANS ${ATLAS_ENABLE_TRANS} )
endif()
if( NOT DEFINED ENABLE_ECTRANS AND DEFINED ENABLE_TRANS )
  ecbuild_warn("Atlas option ENABLE_TRANS is deprecated in favour of ENABLE_ECTRANS")
  set( ENABLE_ECTRANS ${ENABLE_TRANS} )
endif()
if( DEFINED ATLAS_ENABLE_ECTRANS )
  set( ENABLE_ECTRANS ${ATLAS_ENABLE_ECTRANS} )
endif()

set( atlas_HAVE_PACKAGE_ECTRANS 0 )
if( atlas_HAVE_ATLAS_FUNCTIONSPACE AND (ENABLE_ECTRANS OR NOT DEFINED ENABLE_ECTRANS) )
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
        set( atlas_HAVE_PACKAGE_ECTRANS 1 )
    else()
        find_package( transi 0.8 QUIET )
    endif()
endif()
ecbuild_add_option( FEATURE ECTRANS
                    DESCRIPTION "Support for IFS spectral transforms"
                    CONDITION atlas_HAVE_ATLAS_FUNCTIONSPACE AND transi_FOUND )

endif()
