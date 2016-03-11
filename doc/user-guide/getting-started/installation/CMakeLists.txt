ecbuild_add_executable(
  TARGET      atlas_c-hello-world
  SOURCES     hello-world.cc
  LIBS        atlas 
  NOINSTALL )

ecbuild_add_executable(
  TARGET      atlas_f-hello-world
  SOURCES     hello-world.F90
  LIBS        atlas_f 
  NOINSTALL )

set_target_properties(atlas_c-hello-world 
                      atlas_f-hello-world 
                      PROPERTIES EXCLUDE_FROM_ALL true )
