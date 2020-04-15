if( HAVE_TESTS )
  ecbuild_add_option( FEATURE EXPORT_TESTS DEFAULT OFF )
endif()

ecbuild_add_test( TARGET atlas_test_build_hello_world
    COMMAND     ${PROJECT_SOURCE_DIR}/doc/example-projects/build_hello_world.sh
    CONDITION   atlas_HAVE_EXPORT_TESTS
    ENVIRONMENT 
        atlas_DIR=${PROJECT_BINARY_DIR}
        CTEST_ARGS=-VV
        MAKE_ARGS=VERBOSE=1
)

ecbuild_add_test( TARGET atlas_test_build_hello_world_fortran
    COMMAND     ${PROJECT_SOURCE_DIR}/doc/example-projects/build_hello_world_fortran.sh
    CONDITION   atlas_HAVE_EXPORT_TESTS AND atlas_HAVE_FORTRAN
    ENVIRONMENT 
        atlas_DIR=${PROJECT_BINARY_DIR}
        CTEST_ARGS=-VV
        MAKE_ARGS=VERBOSE=1
)
