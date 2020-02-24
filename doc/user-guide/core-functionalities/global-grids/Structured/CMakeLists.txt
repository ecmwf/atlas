ecbuild_add_executable(
  TARGET      atlas_c-global-grids-Structured
  SOURCES     global-grids-Structured.cc
  LIBS        atlas
  NOINSTALL )

set_target_properties(atlas_c-global-grids-Structured
                      PROPERTIES EXCLUDE_FROM_ALL ${ENABLE_TESTS} )

if( atlas_HAVE_FORTRAN )

ecbuild_add_executable(
  TARGET      atlas_f-global-grids-Structured
  SOURCES     global-grids-Structured.F90
  LIBS        atlas_f
  LINKER_LANGUAGE Fortran
  NOINSTALL )

set_target_properties(atlas_f-global-grids-Structured
                      PROPERTIES EXCLUDE_FROM_ALL ${ENABLE_TESTS} )

endif()
