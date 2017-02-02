
if( HAVE_FCTEST )

  add_fctest( TARGET atlas_fctest_griddistribution
    LINKER_LANGUAGE Fortran
    SOURCES         fctest_griddistribution.F90
    LIBS            atlas_f
  )

  add_fctest( TARGET atlas_fctest_state
    LINKER_LANGUAGE Fortran
    SOURCES         fctest_state.F90
    LIBS            atlas_f
  )

endif()

ecbuild_add_test( TARGET atlas_test_state
  BOOST
  SOURCES  test_state.cc
  LIBS     atlas
)


ecbuild_add_test( TARGET atlas_test_grids
  BOOST
  SOURCES   test_grids.cc
  LIBS      atlas
)

ecbuild_add_test(   TARGET      atlas_test_field
    BOOST
    SOURCES     test_field.cc
    LIBS        atlas
)
