ecbuild_add_executable(
  TARGET      atlas_c-global-grids-Unstructured
  SOURCES     global-grids-Unstructured.cc
  LIBS        atlas
  NOINSTALL )

set_target_properties(atlas_c-global-grids-Unstructured
                      PROPERTIES EXCLUDE_FROM_ALL ${ENABLE_TESTS} )

