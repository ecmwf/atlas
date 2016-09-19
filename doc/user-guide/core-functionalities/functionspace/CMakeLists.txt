ecbuild_add_executable(
  TARGET      atlas_c-NodeColumns
  SOURCES     NodeColumns.cc
  LIBS        atlas
  NOINSTALL )

ecbuild_add_executable(
  TARGET      atlas_c-StructuredColumns
  SOURCES     StructuredColumns.cc
  LIBS        atlas
  NOINSTALL )

set_target_properties(atlas_c-NodeColumns
                      atlas_c-StructuredColumns
                      PROPERTIES EXCLUDE_FROM_ALL true )

if( ATLAS_HAVE_FORTRAN )

ecbuild_add_executable(
  TARGET      atlas_f-NodeColumns
  SOURCES     NodeColumns.F90
  LIBS        atlas_f
  NOINSTALL )

set_target_properties(atlas_f-NodeColumns
                      PROPERTIES EXCLUDE_FROM_ALL true )

endif()
