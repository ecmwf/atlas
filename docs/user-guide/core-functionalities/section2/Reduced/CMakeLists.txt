ecbuild_add_executable(
  TARGET      atlas_c-section2-Reduced
  SOURCES     section2-Reduced.cc
  LIBS        atlas
  NOINSTALL )

ecbuild_add_executable(
  TARGET      atlas_f-section2-Reduced
  SOURCES     section2-Reduced.F90
  LIBS        atlas_f
  NOINSTALL )

set_target_properties( atlas_c-section2-Reduced
                       atlas_f-section2-Reduced
                       PROPERTIES EXCLUDE_FROM_ALL true )

