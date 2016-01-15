ecbuild_add_executable(
  TARGET      atlas_c-section2-Gaussian
  SOURCES     section2-Gaussian.cc
  LIBS        atlas
  NOINSTALL )

set_target_properties( atlas_c-section2-Gaussian PROPERTIES EXCLUDE_FROM_ALL true )

