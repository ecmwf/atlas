ecbuild_add_executable(
  TARGET      atlas_c-nodes
  SOURCES     nodes.cc
  LIBS        atlas
  NOINSTALL )

ecbuild_add_executable(
  TARGET      atlas_f-nodes
  SOURCES     nodes.F90
  LIBS        atlas_f
  NOINSTALL )

ecbuild_add_executable(
  TARGET      atlas_c-reducedGridPoint
  SOURCES     reducedGridPoint.cc
  LIBS        atlas
  NOINSTALL )

set_target_properties(atlas_c-nodes
                      atlas_f-nodes
                      atlas_c-reducedGridPoint
                      PROPERTIES EXCLUDE_FROM_ALL true )


