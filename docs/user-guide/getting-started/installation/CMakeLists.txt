ecbuild_add_executable(
  TARGET      atlas_c-usage-example
  SOURCES     usage-example.cc
  LIBS        atlas 
  NOINSTALL )

ecbuild_add_executable(
  TARGET      atlas_f-usage-example
  SOURCES     usage-example.F90
  LIBS        atlas_f 
  NOINSTALL )

set_target_properties(atlas_c-usage-example 
                      atlas_f-usage-example 
                      PROPERTIES EXCLUDE_FROM_ALL true )
