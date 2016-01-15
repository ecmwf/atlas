ecbuild_add_executable(
  TARGET      atlas_c-section1
  SOURCES     section1.cc
  LIBS        atlas 
  NOINSTALL )

ecbuild_add_executable(
  TARGET      atlas_f-section1
  SOURCES     section1.F90
  LIBS        atlas_f 
  NOINSTALL )

set_target_properties( atlas_c-section1 atlas_f-section1 PROPERTIES EXCLUDE_FROM_ALL true )
