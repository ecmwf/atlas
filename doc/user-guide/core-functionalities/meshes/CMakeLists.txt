ecbuild_add_executable(
  TARGET      atlas_c-meshes-Structured
  SOURCES     meshes-Structured.cc
  LIBS        atlas
  NOINSTALL )

set_target_properties(atlas_c-meshes-Structured
                      PROPERTIES EXCLUDE_FROM_ALL ${ENABLE_TESTS} )

if( ATLAS_HAVE_FORTRAN )

ecbuild_add_executable(
  TARGET      atlas_f-meshes-Structured
  SOURCES     meshes-Structured.F90
  LIBS        atlas_f
  NOINSTALL )

set_target_properties(atlas_f-meshes-Structured
                      PROPERTIES EXCLUDE_FROM_ALL ${ENABLE_TESTS} )

endif()
