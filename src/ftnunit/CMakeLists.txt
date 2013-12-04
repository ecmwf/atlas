
add_library( ftnunit
  ftnunit.f90
  table_verification.f90 
)


#ecbuild_add_test( TARGET   test_datastruct
#                  TYPE STATIC
#                  SOURCES  tests/test_datastruct.f90
#                  LIBS     datastruct ftnunit
#                  COMMAND  ${CMAKE_SOURCE_DIR}/src/testing/ftnunit/runtests.sh
#                  ARGS     ${CMAKE_CURRENT_BINARY_DIR}/test_datastruct )
