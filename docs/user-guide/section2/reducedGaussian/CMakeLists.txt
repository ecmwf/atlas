ecbuild_add_executable(
  TARGET      atlas_c-section2-reducedGaussian
  SOURCES     section2-reducedGaussian.cc
  LIBS        atlas
  NOINSTALL )

set_target_properties( atlas_c-section2-reducedGaussian PROPERTIES EXCLUDE_FROM_ALL true )

