ecbuild_add_executable(
  TARGET      atlas_c-global-grids-Reduced
  SOURCES     global-grids-Reduced.cc
  LIBS        atlas
  NOINSTALL )

ecbuild_add_executable(
  TARGET      atlas_f-global-grids-Reduced
  SOURCES     global-grids-Reduced.F90
  LIBS        atlas_f
  NOINSTALL )

set_target_properties(atlas_c-global-grids-Reduced
                      atlas_f-global-grids-Reduced
                      PROPERTIES EXCLUDE_FROM_ALL true )

