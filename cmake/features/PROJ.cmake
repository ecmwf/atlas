
### Proj

ecbuild_find_package( NAME PROJ4 QUIET )
ecbuild_add_option( FEATURE PROJ
                    DESCRIPTION "PROJ-based projections"
                    DEFAULT OFF
                    CONDITION PROJ4_FOUND )
