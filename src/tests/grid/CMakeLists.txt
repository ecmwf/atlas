
if( HAVE_FCTEST )

  add_fctest( TARGET atlas_fctest_griddistribution
    MPI             1
    LINKER_LANGUAGE Fortran
    CONDITION       ON
    SOURCES         fctest_griddistribution.F90
    LIBS            atlas_f
  )

  add_fctest( TARGET atlas_fctest_state
    MPI             1
    LINKER_LANGUAGE Fortran
    SOURCES         fctest_state.F90
    LIBS            atlas_f
  )

endif()

ecbuild_add_test( TARGET atlas_test_state
  BOOST
  MPI      1
  SOURCES  test_state.cc
  LIBS     atlas
)


ecbuild_add_test( TARGET atlas_test_grids
  BOOST
  MPI       1
  SOURCES   test_grids.cc
  LIBS      atlas
)

ecbuild_add_test( TARGET atlas_test_grid_ptr
  BOOST
  SOURCES   test_grid_ptr.cc
  LIBS      atlas
)

ecbuild_add_test(   TARGET      atlas_test_field
    BOOST
    SOURCES     test_field.cc
    LIBS        atlas
    MPI         1 )

