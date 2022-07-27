### Proj

# From proj 9 onwards, it is guaranteed that it was built/installed using CMake and the CMake export is available.
# Then we can safely remove the file FindPROJ.cmake from atlas
# and the following atrocity could be simply replaced with
#
#     ecbuild_add_option( FEATURE PROJ
#                         DESCRIPTION "PROJ-based projections"
#                         DEFAULT OFF
#                         REQUIRED_PACKAGES PROJ )
#

if( ENABLE_PROJ )
    ecbuild_find_package_search_hints( NAME PROJ )

# 1) Try to find PROJ the CMake way (proj >= 7)
    find_package( PROJ CONFIG )
    if( PROJ_FOUND )
        ecbuild_debug("Found PROJ via CONFIG")
    else()
# 2) Try to find PROJ4 the CMake way (proj 6)
        find_package( PROJ4 CONFIG )
        if( PROJ4_FOUND )
            ecbuild_debug("Found PROJ4 via CONFIG")
            set( PROJ_FOUND        ${PROJ4_FOUND} )
            set( PROJ_LIBRARIES    ${PROJ4_LIBRARIES} )
            set( PROJ_INCLUDE_DIRS ${PROJ4_INCLUDE_DIRS} )
            set( PROJ_VERSION      ${PROJ4_VERSION} )
        endif()
    endif()
    if( NOT PROJ_FOUND )
# 3) Try to find PROJ via FindPROJ.cmake provided within atlas (proj < 9)
        find_package( PROJ MODULE )
        if( PROJ_FOUND )
            ecbuild_debug("Found PROJ via FindPROJ.cmake")
        endif()
    endif()
else()
    ecbuild_debug("Skipping search for PROJ as ENABLE_PROJ=OFF (default)")
endif()

if( PROJ_FOUND )
  ecbuild_info( "Found PROJ (${PROJ_VERSION}):")
  ecbuild_info( "    PROJ_LIBRARIES    : ${PROJ_LIBRARIES}")
  ecbuild_info( "    PROJ_INCLUDE_DIRS : ${PROJ_INCLUDE_DIRS}")
endif()

ecbuild_add_option( FEATURE PROJ
                    DESCRIPTION "PROJ-based projections"
                    DEFAULT OFF
                    CONDITION PROJ_FOUND )

if( NOT HAVE_PROJ )
  unset( PROJ_LIBRARIES )
  unset( PROJ_INCLUDE_DIRS )
endif()
