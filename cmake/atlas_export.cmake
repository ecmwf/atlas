################################################################################
# export package info

if( TARGET atlas_f )
  list( APPEND ATLAS_LIBRARIES atlas_f )
endif()
list( APPEND ATLAS_LIBRARIES  atlas )

################################################################################
# pkg-config

ecbuild_add_option( FEATURE PKGCONFIG DESCRIPTION "Atlas pkgconfig" )
set( ATLAS_URL "https://software.ecmwf.int/wiki/display/ATLAS" )
set( ATLAS_DESCRIPTION "Atlas framework for parallel mesh datastructures" )

if( atlas_HAVE_PKGCONFIG )
  ecbuild_pkgconfig(
    NAME             atlas
    LIBRARIES        ${ATLAS_LIBRARIES}
  )

  ecbuild_pkgconfig(
    NAME             atlas-c++
    LANGUAGES        CXX
    LIBRARIES        atlas
  )

  if( atlas_HAVE_FORTRAN )
    ecbuild_pkgconfig(
      NAME             atlas-fortran
      LANGUAGES        Fortran
      LIBRARIES        atlas_f
      NO_PRIVATE_INCLUDE_DIRS
    )
  endif()
endif()

################################################################################
# finalize

ecbuild_add_resources(
    TARGET atlas-others
    SOURCES_PACK
        README.md
        CHANGELOG.md
        LICENSE
)

if( atlas_HAVE_FORTRAN AND ECBUILD_INSTALL_FORTRAN_MODULES )
    install( DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/${CMAKE_CFG_INTDIR}
             DESTINATION module/atlas
             COMPONENT modules )
endif()

set( atlas_REQUIRES_PRIVATE_DEPENDENCIES FALSE )
get_target_property( target_build_type atlas TYPE )
if( target_build_type STREQUAL STATIC_LIBRARY )
  set( atlas_REQUIRES_PRIVATE_DEPENDENCIES TRUE )
endif()

include( atlas_ecbuild2_compatibility )

ecbuild_install_project( NAME Atlas )

