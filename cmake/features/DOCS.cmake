################################################################################
# documentation
if( ENABLE_DOCS )
  find_package(LATEX COMPONENTS PDFLATEX BIBTEX OPTIONAL_COMPONENTS MAKEINDEX HTLATEX)
endif()
ecbuild_add_option( FEATURE DOCS
                    DESCRIPTION "Atlas documentation"
                    DEFAULT OFF
                    CONDITION LATEX_FOUND )
