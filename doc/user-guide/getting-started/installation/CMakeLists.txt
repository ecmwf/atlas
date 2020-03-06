ecbuild_add_executable(
  TARGET      atlas_c-hello-world
  SOURCES     hello-world.cc
  LIBS        atlas
  NOINSTALL )

set_target_properties( atlas_c-hello-world
                       PROPERTIES EXCLUDE_FROM_ALL ${ENABLE_TESTS} )

if( atlas_HAVE_FORTRAN )
  ecbuild_add_executable(
    LINKER_LANGUAGE Fortran
    TARGET      atlas_f-hello-world
    SOURCES     hello-world.F90
    LIBS        atlas_f
    NOINSTALL )
  set_target_properties( atlas_f-hello-world
                         PROPERTIES EXCLUDE_FROM_ALL ${ENABLE_TESTS} )
endif()
