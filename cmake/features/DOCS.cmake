################################################################################
# documentation
if( ENABLE_DOCS ) 
  find_package(Latex)
endif()
ecbuild_add_option( FEATURE DOCS
                    DESCRIPTION "Atlas documentation"
                    DEFAULT OFF
                    CONDITION Latex_FOUND )
