

ecbuild_add_test( TARGET atlas_test_repartition
  SOURCES   test_repartition.cc
  MPI       8
  LIBS      atlas
  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
)
