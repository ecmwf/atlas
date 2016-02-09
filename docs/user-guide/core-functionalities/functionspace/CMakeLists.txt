ecbuild_add_executable(
  TARGET      atlas_c-functionspace
  SOURCES     functionspace.cc
  LIBS        atlas
  NOINSTALL )

ecbuild_add_executable(
  TARGET      atlas_f-functionspace
  SOURCES     functionspace.F90
  LIBS        atlas_f
  NOINSTALL )

set_target_properties(atlas_c-functionspace
                      atlas_f-functionspace
                      PROPERTIES EXCLUDE_FROM_ALL true )


