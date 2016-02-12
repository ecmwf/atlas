add_subdirectory( user-guide )

add_custom_target( docs )
add_dependencies ( docs 
                   atlas-user-guide-pdf
                   atlas-user-guide-html
                   atlas_c-usage-example
                   atlas_f-usage-example
                   atlas_c-global-grids-Reduced
                   atlas_f-global-grids-Reduced
                   atlas_c-meshes-Reduced
                   atlas_f-meshes-Reduced
                   atlas_c-fields
                   atlas_f-fields
                   atlas_c-fields-on-grid
                   atlas_f-fields-on-grid
                   atlas_c-nodes
                   atlas_f-nodes )
