add_subdirectory( section1 )
add_subdirectory( section2 )

set( USERGUIDESRC ${CMAKE_CURRENT_SOURCE_DIR} )
set( USERGUIDE ${CMAKE_BINARY_DIR}/atlas/docs/user-guide )


# HTML GENERATION
file( MAKE_DIRECTORY ${USERGUIDE}/html )

add_custom_target( atlas-user-guide-html
    export TEXINPUTS=${CMAKE_SOURCE_DIR}//:${USERGUIDESRC}//: &&
    ${HTLATEX} ${USERGUIDESRC}/user-guide.tex 
    "${USERGUIDESRC}/styling.cfg,html,3,next,NoFonts"
    WORKING_DIRECTORY ${USERGUIDE}/html )

# If tex4ht successful, create img dir and copy images across
FILE(GLOB_RECURSE imgfiles "*/img/*.png" "*/img/*.jpg")
ADD_CUSTOM_COMMAND(TARGET user-guide-html
    POST_BUILD COMMAND ${CMAKE_COMMAND} -E make_directory ${USERGUIDE}/html/img)
FOREACH(img ${imgfiles})
    ADD_CUSTOM_COMMAND(TARGET user-guide-html
        POST_BUILD COMMAND 
            ${CMAKE_COMMAND} -E copy ${img} ${USERGUIDE}/html/img)
ENDFOREACH()

FILE(GLOB_RECURSE pdffiles "*/img/*.pdf")
FIND_PROGRAM(CONVERT convert)
FOREACH(pdf ${pdffiles})
    GET_FILENAME_COMPONENT(BASENAME ${pdf} NAME_WE)
    ADD_CUSTOM_COMMAND(TARGET user-guide-html
        POST_BUILD COMMAND
        ${CONVERT} ${pdf} ${USERGUIDE}/html/img/${BASENAME}.png)
ENDFOREACH()

# PDF GENERATION
add_custom_target( atlas-user-guide-pdf
    export TEXINPUTS=${CMAKE_SOURCE_DIR}//: && 
    ${PDFLATEX} --output-directory ${USERGUIDE} ${USERGUIDESRC}/user-guide.tex
    COMMAND export TEXINPUTS=${CMAKE_SOURCE_DIR}//: && 
    ${PDFLATEX} --output-directory ${USERGUIDE} ${USERGUIDESRC}/user-guide.tex
    WORKING_DIRECTORY ${USERGUIDESRC} )
