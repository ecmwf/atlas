# (C) Copyright 2013 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

if(atlas_HAVE_DOCS)

if(atlas_HAVE_FORTRAN)
  add_subdirectory(example-fortran)
endif()

add_subdirectory(user-guide)

list( APPEND DOC_CODE_TARGETS
  atlas_c-hello-world
  atlas_c-global-grids-Structured
  atlas_c-global-grids-Unstructured
  atlas_c-meshes-Structured
  atlas_c-fields
  atlas_c-fields-on-grid
  atlas_c-NodeColumns
  atlas_c-StructuredColumns
)
if( atlas_HAVE_FORTRAN )
  list( APPEND DOC_CODE_TARGETS
    atlas_f-hello-world
    atlas_f-global-grids-Structured
    atlas_f-meshes-Structured
    atlas_f-fields
    atlas_f-fields-on-grid
    atlas_f-NodeColumns
  )
endif()

add_custom_target(atlas_doc_code)
add_dependencies (atlas_doc_code
                  ${DOC_CODE_TARGETS} )

add_custom_target(atlas_doc_pdf)
add_dependencies (atlas_doc_pdf
                  atlas-user-guide-pdf )

if( LATEX_HTLATEX_FOUND )
  add_custom_target(atlas_doc_html)
  add_dependencies (atlas_doc_html
                    atlas-user-guide-html )
endif()

add_custom_target(atlas_doc)
add_dependencies (atlas_doc
                  atlas_doc_code
                  atlas_doc_pdf )

if( PROJECT_NAME STREQUAL CMAKE_PROJECT_NAME )

  add_custom_target(doc)
  add_dependencies(doc atlas_doc)

endif()


if( NOT DEFINED ATLAS_DOC_VERSION )
  set( ATLAS_DOC_VERSION ${ATLAS_VERSION_STR} )
endif()

if( NOT DEFINED ATLAS_DOXYGEN_GENERATOR )
  set( ATLAS_DOXYGEN_GENERATOR "stock")
endif()

if( ATLAS_DOXYGEN_GENERATOR STREQUAL "m.css" )
  set( ATLAS_DOXYFILE Doxyfile-mcss )
  set( ATLAS_DOXYGEN_EXECUTABLE doxygen.py )
else()
  set( ATLAS_DOXYFILE Doxyfile-stock )
endif()

if( NOT DEFINED ATLAS_DOXYGEN_EXECUTABLE )
  find_package(Doxygen COMPONENTS dot)
  if( DOXYGEN_FOUND )
    set( ATLAS_DOXYGEN_EXECUTABLE ${DOXYGEN_EXECUTABLE} )
  else()
    set( ATLAS_DOXYGEN_EXECUTABLE doxygen )
  endif()
endif()

foreach( doxyfile Doxyfile-default Doxyfile-custom Doxyfile-mcss Doxyfile-stock )
  execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${CMAKE_CURRENT_SOURCE_DIR}/${doxyfile}
      ${CMAKE_CURRENT_BINARY_DIR}/${doxyfile} )
endforeach()

configure_file(Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)

add_custom_target(
  atlas_doxygen
  COMMAND ${ATLAS_DOXYGEN_EXECUTABLE} Doxyfile
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  COMMENT "Generating Doxygen documentation" VERBATIM
)

set_property(TARGET atlas_doxygen PROPERTY EXCLUDE_FROM_ALL TRUE)

endif()
