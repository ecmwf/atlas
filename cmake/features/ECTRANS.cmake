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

macro( make_transi_alias target )
  if( NOT TARGET transi )
    if( CMAKE_VERSION VERSION_LESS 3.18 )
      # Before CMake 3.18 it is not possible to alias a non-global imported target
      # Make the import global. Warning, this may break further find_package
      get_target_property( ${target}_IMPORTED ${target} IMPORTED )
      if( transi_dp_IMPORTED )
        set_target_properties( ${target} PROPERTIES IMPORTED_GLOBAL TRUE)
      endif()
    endif()
    add_library( transi ALIAS ${target} )
  endif()
endmacro()

set( atlas_HAVE_PACKAGE_ECTRANS 0 )
if( atlas_HAVE_ATLAS_FUNCTIONSPACE AND (ENABLE_ECTRANS OR NOT DEFINED ENABLE_ECTRANS) )
    find_package( ectrans 1.1 COMPONENTS transi double QUIET )
    if( TARGET transi_dp OR TARGET transi_gpu_dp )
        set( transi_FOUND TRUE )
        set( atlas_HAVE_PACKAGE_ECTRANS 1 )
    else()
        find_package( transi 0.8 QUIET )
    endif()
endif()
ecbuild_add_option( FEATURE ECTRANS
                    DESCRIPTION "Support for ectrans spectral transforms"
                    CONDITION atlas_HAVE_ATLAS_FUNCTIONSPACE AND transi_FOUND )
ecbuild_add_option( FEATURE ECTRANS_GPU DEFAULT OFF
                    DESCRIPTION "Support for ectrans spectral transforms"
                    CONDITION HAVE_ECTRANS AND TARGET transi_gpu_dp )
if( HAVE_ECTRANS_GPU AND TARGET transi_gpu_dp)
  ecbuild_warn("ectrans GPU version")
  make_transi_alias( transi_gpu_dp )
elseif( HAVE_ECTRANS AND TARGET transi_dp )
  ecbuild_warn("ectrans CPU version")
  make_transi_alias( transi_dp )
endif()

endif()
