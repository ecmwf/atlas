ecbuild_add_executable(
  TARGET      atlas_c-meshes-Reduced
  SOURCES     meshes-Reduced.cc
  LIBS        atlas
  NOINSTALL )

ecbuild_add_executable(
  TARGET      atlas_f-meshes-Reduced
  SOURCES     meshes-Reduced.F90
  LIBS        atlas_f
  NOINSTALL )

set_target_properties(atlas_c-meshes-Reduced 
                      atlas_f-meshes-Reduced
                      PROPERTIES EXCLUDE_FROM_ALL true )


