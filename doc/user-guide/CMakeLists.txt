# (C) Copyright 2013 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

configure_file( version.tex.in version.tex )

add_subdirectory(getting-started)
add_subdirectory(core-functionalities)

set(USERGUIDESRC ${CMAKE_CURRENT_SOURCE_DIR})
set(USERGUIDE ${CMAKE_CURRENT_BINARY_DIR})



# HTML GENERATION =============================================================
if( LATEX_HTLATEX_FOUND )
file(MAKE_DIRECTORY ${USERGUIDE}/html)

  add_custom_target(atlas-user-guide-html
    export TEXINPUTS=${CMAKE_SOURCE_DIR}//:${USERGUIDESRC}//:${CMAKE_CURRENT_BINARY_DIR}//: &&
    ${HTLATEX_COMPILER} ${USERGUIDESRC}/user-guide.tex
    "${USERGUIDESRC}/styling.cfg,html,3,next,NoFonts"
    WORKING_DIRECTORY ${USERGUIDE}/html )

# If tex4ht successful, create img dir and copy images across
# .png and .jpg images --------------------------------------------------------
  FILE(GLOB_RECURSE imgfiles
      "introduction/imgs/*.png" "introduction/imgs/*.jpg"
      "getting-started/*/imgs/*.png" "getting-started/*/imgs/*jpg"
      "core-functionalities/*/imgs/*.png" "core-functionalities/*/imgs/*jpg"
      "applications/*/imgs/*.png" "applications/*/imgs/*jpg")

  add_custom_command(
     TARGET atlas-user-guide-html
     POST_BUILD COMMAND ${CMAKE_COMMAND} -E make_directory ${USERGUIDE}/html/imgs)

  FOREACH(img ${imgfiles})
     ADD_CUSTOM_COMMAND(TARGET atlas-user-guide-html
         POST_BUILD COMMAND
         ${CMAKE_COMMAND} -E copy ${img} ${USERGUIDE}/html/imgs)
  ENDFOREACH()
# -----------------------------------------------------------------------------

# .pdf images -----------------------------------------------------------------
  FILE(GLOB_RECURSE pdffiles
      "introduction/imgs/*.pdf"
      "getting-started/*/imgs/*.pdf"
      "core-functionalities/*/imgs/*.pdf"
      "applications/*/imgs/*.pdf")

  FIND_PROGRAM(CONVERT convert)

  FOREACH(pdf ${pdffiles})
     GET_FILENAME_COMPONENT(BASENAME ${pdf} NAME_WE)
     ADD_CUSTOM_COMMAND(
         TARGET atlas-user-guide-html
         POST_BUILD COMMAND
         ${CONVERT} ${pdf} ${USERGUIDE}/html/imgs/${BASENAME}.png)
  ENDFOREACH()
endif( LATEX_HTLATEX_FOUND )

# =============================================================================



# PDF GENERATION ==============================================================
add_custom_target(atlas-user-guide-pdf
   export TEXINPUTS=${CMAKE_SOURCE_DIR}//:${CMAKE_CURRENT_BINARY_DIR}//: &&
   ${PDFLATEX_COMPILER} --output-directory ${USERGUIDE} ${USERGUIDESRC}/user-guide.tex
   COMMAND export TEXINPUTS=${CMAKE_SOURCE_DIR}//:${CMAKE_CURRENT_BINARY_DIR}//: &&
   ${PDFLATEX_COMPILER} --output-directory ${USERGUIDE} ${USERGUIDESRC}/user-guide.tex
   WORKING_DIRECTORY ${USERGUIDESRC} )
# =============================================================================
