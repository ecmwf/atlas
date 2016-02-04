add_subdirectory( user-guide )

add_custom_target( docs )
add_dependencies ( docs 
                   atlas-user-guide-pdf
                   atlas-user-guide-html
                   atlas_c-section1
                   atlas_f-section1 
                   atlas_c-section2-Reduced
                   atlas_f-section2-Reduced )
