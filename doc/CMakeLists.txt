# (C) Copyright 2013 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

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
if( ATLAS_HAVE_FORTRAN )
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

add_custom_target(atlas_doc_html)
add_dependencies (atlas_doc_html
                  atlas-user-guide-html )

add_custom_target(atlas_doc)
add_dependencies (atlas_doc
                  atlas_doc_code
                  atlas_doc_pdf )

if( PROJECT_NAME STREQUAL CMAKE_PROJECT_NAME )

  add_custom_target(doc)
  add_dependencies(doc atlas_doc)

endif()


