################################################################################
# documentation
if( ENABLE_DOCS )
  find_package(LATEX REQUIRED COMPONENTS PDFLATEX BIBTEX MAKEINDEX HTLATEX)
endif()
ecbuild_add_option( FEATURE DOCS
                    DESCRIPTION "Atlas documentation"
                    DEFAULT OFF
                    CONDITION LATEX_FOUND )
