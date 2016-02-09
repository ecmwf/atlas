ecbuild_add_executable(
  TARGET      atlas_c-fields
  SOURCES     fields.cc
  LIBS        atlas
  NOINSTALL )

ecbuild_add_executable(
  TARGET      atlas_c-fields-on-grid
  SOURCES     fields-on-grid.cc
  LIBS        atlas
  NOINSTALL )

ecbuild_add_executable(
  TARGET      atlas_f-fields
  SOURCES     fields.F90
  LIBS        atlas_f
  NOINSTALL )

ecbuild_add_executable(
  TARGET      atlas_f-fields-on-grid
  SOURCES     fields-on-grid.F90
  LIBS        atlas_f
  NOINSTALL )

set_target_properties(atlas_c-fields
                      atlas_f-fields
                      atlas_c-fields-on-grid
                      atlas_f-fields-on-grid
                      PROPERTIES EXCLUDE_FROM_ALL true )


