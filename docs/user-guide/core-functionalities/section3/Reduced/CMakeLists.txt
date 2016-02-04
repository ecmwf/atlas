ecbuild_add_executable(
  TARGET      atlas_c-section3-Reduced
  SOURCES     section3-Reduced.cc
  LIBS        atlas
  NOINSTALL )

set_target_properties( atlas_c-section3-Reduced PROPERTIES EXCLUDE_FROM_ALL true )

