ecbuild_add_test( TARGET test_datastruct
  LINKER_LANGUAGE Fortran
  SOURCES  test_datastruct.F90
  LIBS     ftnunit datastruct
  COMMAND  ${CMAKE_SOURCE_DIR}/src/ftnunit/runtests.sh
  ARGS     ${CMAKE_CURRENT_BINARY_DIR}/test_datastruct
  CONDITION  OFF
)

ecbuild_add_test( TARGET test_gmsh
  SOURCES  test_gmsh.cpp
  LIBS     atlas
)
