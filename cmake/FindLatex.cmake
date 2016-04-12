# © Copyright 1996-2016 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

###############################################################################

# find latex libraries
FIND_PROGRAM (PDFLATEX  
              NAMES pdflatex 
              PATHS /usr/local/share/apps/TeXLive/2014/bin/x86_64-linux/)
FIND_PROGRAM (BIBTEX  
              NAMES bibtex 
              PATHS /usr/local/share/apps/TeXLive/2014/bin/x86_64-linux/)
FIND_PROGRAM (MAKEINDEX  
              NAMES makeindex 
              PATHS /usr/local/share/apps/TeXLive/2014/bin/x86_64-linux/)
FIND_PROGRAM (HTLATEX  
              NAMES htlatex 
              PATHS /usr/local/share/apps/TeXLive/2014/bin/x86_64-linux/)

if ( PDFLATEX AND BIBTEX AND MAKEINDEX AND HTLATEX ) 
    set( LATEX_FOUND TRUE )
endif()


mark_as_advanced(PDFLATEX BIBTEX MAKEINDEX HTLATEX)


