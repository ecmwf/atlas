add_subdirectory(user-guide)

list( APPEND DOC_CODE_TARGETS 
  atlas_c-hello-world
  atlas_c-global-grids-Reduced
  atlas_c-meshes-Reduced
  atlas_c-fields
  atlas_c-fields-on-grid
  atlas_c-nodes
  atlas_c-reducedGridPoint
)
if( ATLAS_HAVE_FORTRAN )
  list( APPEND DOC_CODE_TARGETS 
    atlas_f-hello-world
    atlas_f-global-grids-Reduced
    atlas_f-meshes-Reduced
    atlas_f-fields
    atlas_f-fields-on-grid
    atlas_f-nodes
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


